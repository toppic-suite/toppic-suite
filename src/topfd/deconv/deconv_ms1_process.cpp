//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//

#include "topfd/deconv/deconv_ms1_process.hpp"

#include <cstddef>
#include <string>

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

//add a namespace to avoid duplicated method names
namespace deconv_ms1_process {

void deconvMsOne(const MzmlMsGroupPtr &ms_group_ptr, 
                 const TopfdParaPtr &topfd_para_ptr,  
                 const MsAlignWriterPtrVec &ms1_writer_ptr_vec,
                 const SimpleThreadPoolPtr &pool_ptr,
                 const MzmlMsSqlWriterPtr &sql_writer_ptr) {
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
  MatchEnvPtrVec prec_envs = deconv_prec_win::deconvPrecWinForMsGroup(ms_group_ptr, 
                                                                      topfd_para_ptr->getMaxMass(),
                                                                      topfd_para_ptr->getMaxCharge(),
                                                                      base_inte, min_ref_inte); 

  // Obtain EnvCNN Score for envelopes
  onnx_env_cnn::computeEnvScores(peak_list, prec_envs); 

  //remove precursor peaks
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
    DeconvSingleSpPtr deconv_ptr 
      = std::make_shared<DeconvSingleSp>(topfd_para_ptr, peak_list, ms_level,
                                         max_mass, max_charge);
    deconv_envs = deconv_ptr->deconv();
  }
  // 4. Merge precursor envelopes and deconvolution envelopes
  MatchEnvPtrVec result_envs;
  result_envs.insert(result_envs.end(), prec_envs.begin(), prec_envs.end());
  result_envs.insert(result_envs.end(), deconv_envs.begin(), deconv_envs.end());
  LOG_DEBUG("result num " << result_envs.size());
  
  // 5. Write to msalign file
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  DeconvMsPtr deconv_ms_ptr = match_env_util::getDeconvMsPtr(header_ptr,
                                                             result_envs);
  

  std::thread::id thread_id = std::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);
  ms1_writer_ptr_vec[writer_id]->writeMs(deconv_ms_ptr);
  
  // 6. write the deconvoluted spectrum to the SQLite database (if enabled)
  if (sql_writer_ptr != nullptr) {
    sql_writer_ptr->writeMs1(ms_ptr, result_envs, base_inte, min_ref_inte);
  }

}

std::function<void()> geneTask(const MzmlMsGroupPtr &ms_group_ptr,
                               const TopfdParaPtr &topfd_para_ptr,
                               const MsAlignWriterPtrVec &ms1_writer_ptr_vec,
                               const SimpleThreadPoolPtr &pool_ptr,
                               const MzmlMsSqlWriterPtr &sql_writer_ptr) {
  return [ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr, sql_writer_ptr]() {
    deconvMsOne(ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr, sql_writer_ptr);
  };
}

} // namespace deconv_ms1_process end

DeconvMs1Process::DeconvMs1Process(const TopfdParaPtr &topfd_para_ptr) {
  topfd_para_ptr_ = topfd_para_ptr;
}

void DeconvMs1Process::process() {
  MzmlMsGroupReaderPtr reader_ptr = 
    std::make_shared<MzmlMsGroupReader>(topfd_para_ptr_->getMzmlFileName(), 
                                        topfd_para_ptr_->getPrecWindowWidth(),
                                        topfd_para_ptr_->getActivation(),
                                        topfd_para_ptr_->getFracId(),
                                        topfd_para_ptr_->isFaims(), 
                                        topfd_para_ptr_->getFaimsVoltage(), 
                                        topfd_para_ptr_->isMissingLevelOne());

  MzmlMsGroupPtr ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  if (ms_group_ptr == nullptr) {
    LOG_ERROR("No spectrum to read in mzML file!");
    return;
  }
  // One SQLite writer shared across the worker threads (it is internally
  // synchronized and batches inserts); created only when SQLite output is on.
  MzmlMsSqlWriterPtr sql_writer_ptr = topfd_para_ptr_->isGeneSql()
      ? std::make_shared<MzmlMsSqlWriter>(topfd_para_ptr_->getSqlDb())
      : nullptr;
  // init thread pool
  int thread_num = topfd_para_ptr_->getThreadNum();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);  
  // init msalign writer vector for multiple threads
  std::string output_base_name = topfd_para_ptr_->getOutputBaseName();
  std::string ms1_msalign_name = output_base_name + "_ms1.msalign";
  MsAlignWriterPtrVec ms1_writer_ptr_vec;
  for (int i = 0; i < thread_num; i++) { 
    MsAlignWriterPtr ms1_ptr 
        = std::make_shared<MsAlignWriter>(ms1_msalign_name + "_" + std::to_string(i));
    ms1_writer_ptr_vec.push_back(ms1_ptr);
  }
  // counter for processed spectra
  int spec_cnt = 0;
  // total spectrum number
  int total_spec_num = topfd_para_ptr_->getMs1ScanNum(); 
  while (ms_group_ptr != nullptr) {
    while(pool_ptr->getQueueSize() >= thread_num * 2){
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    pool_ptr->enqueue(deconv_ms1_process::geneTask(ms_group_ptr, topfd_para_ptr_, ms1_writer_ptr_vec, pool_ptr, sql_writer_ptr));
    spec_cnt++;
    std::string msg = deconv_util::updateMsOneMsg(ms_group_ptr->getMsOnePtr()->getMsHeaderPtr(),
                                                  spec_cnt, total_spec_num);
    std::cout << "\r" << msg << std::flush;
    ms_group_ptr = reader_ptr->getNextMsGroupPtr();    
  }
  pool_ptr->shutDown();
  if (sql_writer_ptr != nullptr) {
    sql_writer_ptr->flush();
  }
  for (int i = 0; i < thread_num; i++) {
    ms1_writer_ptr_vec[i] = nullptr;
  }
  deconv_util::mergeMs1MsalignFiles(topfd_para_ptr_, output_base_name);
}

}; // namespace toppic
