// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "topfd/deconv/deconv_ms1_process.hpp"

#include <atomic>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/logger.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/mzml/mzml_ms_sql_writer.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/spec/msalign_thread_merge.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "topfd/deconv/deconv_prec_win.hpp"
#include "topfd/deconv/deconv_single_sp.hpp"
#include "topfd/deconv/deconv_util.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

// add a namespace to avoid duplicated method names
namespace deconv_ms1_process {

// One deconvoluted MS1 spectrum's data, buffered in memory so the slow SQLite
// writes can be deferred and done sequentially after the parallel
// deconvolution finishes, instead of every worker thread contending on the
// writer lock and the SSD.
struct Ms1SqlRecord {
  MzmlMsPtr ms_ptr;
  MatchEnvPtrVec result_envs;
  double base_inte;
  double min_ref_inte;
};

// One record buffer per worker thread (indexed by writer id) so threads append
// without locking, plus an approximate byte counter shared across threads so
// the main thread can cap the buffer's memory use and flush it when full.
struct Ms1SqlBuffer {
  explicit Ms1SqlBuffer(int thread_num) : thread_records(thread_num) {}
  std::vector<std::vector<Ms1SqlRecord>> thread_records;
  std::atomic<std::size_t> bytes{0};
};
using Ms1SqlBufferPtr = std::shared_ptr<Ms1SqlBuffer>;

// Approximate bytes held per buffered peak (the Peak/EnvPeak object, its
// shared_ptr control block, and the vector slot). Used only to bound the
// buffer's memory, so a rough value is sufficient.
constexpr std::size_t BYTES_PER_PEAK = 64;

// Hard ceiling for the in-memory SQL buffer.
constexpr std::size_t MAX_SQL_BUFFER_BYTES =
    4ULL * 1024 * 1024 * 1024;  // 4 GiB

// Start draining at 90% of the ceiling, leaving headroom for the records that
// in-flight worker tasks keep appending while the buffer is being drained, so
// the actual peak stays under MAX_SQL_BUFFER_BYTES.
constexpr std::size_t SQL_BUFFER_FLUSH_BYTES = MAX_SQL_BUFFER_BYTES * 9 / 10;

// Rough memory estimate for one buffered record: the spectrum peaks plus the
// peaks of all its envelopes.
std::size_t estimateRecordBytes(const Ms1SqlRecord& record) {
  std::size_t peak_num = record.ms_ptr->size();
  for (const auto& env : record.result_envs) {
    peak_num += env->getExpEnvPtr()->getPeakNum();
  }
  return peak_num * BYTES_PER_PEAK;
}

void deconvMsOne(const MzmlMsGroupPtr& ms_group_ptr,
                 const TopfdParaPtr& topfd_para_ptr,
                 const MsAlignWriterPtrVec& ms1_writer_ptr_vec,
                 const SimpleThreadPoolPtr& pool_ptr,
                 const Ms1SqlBufferPtr& sql_buffer_ptr) {
  // 1. Store peak intensity
  MzmlMsPtr ms_ptr = ms_group_ptr->getMsOnePtr();
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  std::vector<double> intensities;
  for (std::size_t i = 0; i < peak_list.size(); i++) {
    intensities.push_back(peak_list[i]->getIntensity());
  }
  double base_inte = baseline_util::getBaseLine(intensities);
  double min_ref_inte = base_inte * topfd_para_ptr->getMsOneSnRatio();

  // 2. Deconv envelopes in precursor windows and remove them
  MatchEnvPtrVec prec_envs = deconv_prec_win::deconvPrecWinForMsGroup(
      ms_group_ptr, topfd_para_ptr->getMaxMass(),
      topfd_para_ptr->getMaxCharge(), base_inte, min_ref_inte);

  // Obtain EnvCNN Score for envelopes
  onnx_env_cnn::computeEnvScores(peak_list, prec_envs);

  // remove precursor peaks
  for (std::size_t i = 0; i < prec_envs.size(); i++) {
    ExpEnvPtr env_ptr = prec_envs[i]->getExpEnvPtr();
    for (int p = 0; p < env_ptr->getPeakNum(); p++) {
      if (env_ptr->isExist(p)) {
        int idx = env_ptr->getPeakIdx(p);
        peak_list[idx]->setIntensity(0);
      }
    }
  }
  // 3. Deconv the whole spectrum with filtering
  // get base intensity and min_ref_intensity for sql writing
  MatchEnvPtrVec deconv_envs;
  if (peak_list.size() > 0) {
    int ms_level = 1;
    double max_mass = topfd_para_ptr->getMaxMass();
    int max_charge = topfd_para_ptr->getMaxCharge();
    DeconvSingleSpPtr deconv_ptr = std::make_shared<DeconvSingleSp>(
        topfd_para_ptr, peak_list, ms_level, max_mass, max_charge);
    deconv_envs = deconv_ptr->deconv();
  }
  // Restore the intensities of the removed precursor peaks (zeroed above so
  // the whole-spectrum deconvolution skips them). The peaks are shared with
  // ms_ptr, which is buffered for SQLite output below, so the spectrum must be
  // written with its original intensities.
  for (std::size_t i = 0; i < prec_envs.size(); i++) {
    ExpEnvPtr env_ptr = prec_envs[i]->getExpEnvPtr();
    for (int p = 0; p < env_ptr->getPeakNum(); p++) {
      if (env_ptr->isExist(p)) {
        int idx = env_ptr->getPeakIdx(p);
        peak_list[idx]->setIntensity(intensities[idx]);
      }
    }
  }
  // 4. Merge precursor envelopes and deconvolution envelopes
  MatchEnvPtrVec result_envs;
  result_envs.insert(result_envs.end(), prec_envs.begin(), prec_envs.end());
  result_envs.insert(result_envs.end(), deconv_envs.begin(), deconv_envs.end());
  LOG_DEBUG("result num " << result_envs.size());

  // 5. Write to msalign file
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  DeconvMsPtr deconv_ms_ptr =
      match_env_util::getDeconvMsPtr(header_ptr, result_envs);

  std::thread::id thread_id = std::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);
  ms1_writer_ptr_vec[writer_id]->writeMs(deconv_ms_ptr);

  // 6. Buffer the deconvoluted spectrum in memory for SQLite output (if
  // enabled). Each worker thread appends to its own buffer (indexed by writer
  // id), so there is no lock contention; the records are written to the
  // database sequentially after the thread pool finishes.
  if (sql_buffer_ptr != nullptr) {
    Ms1SqlRecord record{ms_ptr, std::move(result_envs), base_inte,
                        min_ref_inte};
    sql_buffer_ptr->bytes += estimateRecordBytes(record);
    sql_buffer_ptr->thread_records[writer_id].push_back(std::move(record));
  }
}

std::function<void()> geneTask(const MzmlMsGroupPtr& ms_group_ptr,
                               const TopfdParaPtr& topfd_para_ptr,
                               const MsAlignWriterPtrVec& ms1_writer_ptr_vec,
                               const SimpleThreadPoolPtr& pool_ptr,
                               const Ms1SqlBufferPtr& sql_buffer_ptr) {
  return [ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr,
          sql_buffer_ptr]() {
    deconvMsOne(ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr,
                sql_buffer_ptr);
  };
}

}  // namespace deconv_ms1_process

DeconvMs1Process::DeconvMs1Process(const TopfdParaPtr& topfd_para_ptr) {
  topfd_para_ptr_ = topfd_para_ptr;
}

void DeconvMs1Process::process() {
  MzmlMsGroupReaderPtr reader_ptr = std::make_shared<MzmlMsGroupReader>(
      topfd_para_ptr_->getMzmlFileName(), topfd_para_ptr_->getPrecWindowWidth(),
      topfd_para_ptr_->getActivation(), topfd_para_ptr_->getFracId(),
      topfd_para_ptr_->isFaims(), topfd_para_ptr_->getFaimsVoltage(),
      topfd_para_ptr_->isMissingLevelOne());

  MzmlMsGroupPtr ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  if (ms_group_ptr == nullptr) {
    LOG_ERROR("No spectrum to read in mzML file!");
    return;
  }
  // One SQLite writer shared across the worker threads (it is internally
  // synchronized and batches inserts); created only when SQLite output is on.
  MzmlMsSqlWriterPtr sql_writer_ptr =
      topfd_para_ptr_->isGeneSql()
          ? std::make_shared<MzmlMsSqlWriter>(topfd_para_ptr_->getSqlDb())
          : nullptr;
  // init thread pool
  int thread_num = topfd_para_ptr_->getThreadNum();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);
  // Per-thread in-memory buffers for the SQLite records. Writing to SQLite from
  // every worker thread serializes on the writer lock and the SSD; instead each
  // thread buffers its records here (one sub-vector per thread, no locking) and
  // they are written to the database sequentially once deconvolution is done.
  deconv_ms1_process::Ms1SqlBufferPtr sql_buffer_ptr =
      sql_writer_ptr != nullptr
          ? std::make_shared<deconv_ms1_process::Ms1SqlBuffer>(thread_num)
          : nullptr;

  // Wait until the pool has drained: no queued tasks and every worker idle, so
  // no thread is touching the SQL buffer and the main thread can flush it.
  auto wait_until_pool_idle = [&]() {
    while (pool_ptr->getQueueSize() > 0 ||
           pool_ptr->getIdleThreadNum() < thread_num) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  };
  // Write the buffered records to SQLite sequentially and empty the buffer.
  // Must be called only while the pool is idle (see wait_until_pool_idle).
  auto flush_sql_buffer = [&]() {
    for (auto& thread_buffer : sql_buffer_ptr->thread_records) {
      for (const auto& record : thread_buffer) {
        sql_writer_ptr->writeMs1(record.ms_ptr, record.result_envs,
                                 record.base_inte, record.min_ref_inte);
      }
      thread_buffer.clear();
    }
    sql_buffer_ptr->bytes = 0;
  };

  // init msalign writer vector for multiple threads
  std::string output_base_name = topfd_para_ptr_->getOutputBaseName();
  std::string ms1_msalign_name = output_base_name + "_ms1.msalign";
  MsAlignWriterPtrVec ms1_writer_ptr_vec;
  for (int i = 0; i < thread_num; i++) {
    MsAlignWriterPtr ms1_ptr = std::make_shared<MsAlignWriter>(
        ms1_msalign_name + "_" + std::to_string(i));
    ms1_writer_ptr_vec.push_back(ms1_ptr);
  }
  // counter for processed spectra
  int spec_cnt = 0;
  // total spectrum number
  int total_spec_num = topfd_para_ptr_->getMs1ScanNum();
  while (ms_group_ptr != nullptr) {
    while (pool_ptr->getQueueSize() >= thread_num * 2) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    pool_ptr->enqueue(deconv_ms1_process::geneTask(
        ms_group_ptr, topfd_para_ptr_, ms1_writer_ptr_vec, pool_ptr,
        sql_buffer_ptr));
    spec_cnt++;
    std::string msg = deconv_util::updateMsOneMsg(
        ms_group_ptr->getMsOnePtr()->getMsHeaderPtr(), spec_cnt,
        total_spec_num);
    std::cout << "\r" << msg << std::flush;

    // Backpressure: if the in-memory SQL buffer is full, pause deconvolution,
    // drain the buffer to the database, and then resume.
    if (sql_writer_ptr != nullptr &&
        sql_buffer_ptr->bytes.load() >=
            deconv_ms1_process::SQL_BUFFER_FLUSH_BYTES) {
      wait_until_pool_idle();
      flush_sql_buffer();
    }

    ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  }
  pool_ptr->shutDown();
  // Write any remaining buffered records to SQLite sequentially (single-
  // threaded, so no lock contention) and commit the final batch.
  if (sql_writer_ptr != nullptr) {
    flush_sql_buffer();
    sql_writer_ptr->flush();
  }
  for (int i = 0; i < thread_num; i++) {
    ms1_writer_ptr_vec[i] = nullptr;
  }
  deconv_util::mergeMs1MsalignFiles(topfd_para_ptr_, output_base_name);
}

};  // namespace toppic
