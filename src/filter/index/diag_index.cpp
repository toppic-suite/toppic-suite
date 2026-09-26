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

#include "filter/index/diag_index.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <system_error>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace diag_index {

// serialization mutex.
namespace {
std::mutex serial_mutex;

std::string indexDir(const DiagFilterMngPtr& mng_ptr) {
  return mng_ptr->prsm_para_ptr_->getOriDbName() + "_idx" +
         file_util::getFileSeparator();
}
}  // namespace

// The multi_ptm_index of a block is the PRM diagonal MassMatch of its
// proteoforms at mng_ptr's column limits - the same object the one-shift
// stage builds and writes as one_ptm_diag_index when the two managers share
// max_proteoform_mass_ and filter_scale_ (both 50000 and 100 by default), so
// the two files are byte-identical and the second build is redundant.
bool canReuseOnePtmDiagIndex(const DiagFilterMngPtr& mng_ptr,
                             const OnePtmFilterMngPtr& one_ptm_mng_ptr) {
  return one_ptm_mng_ptr != nullptr &&
         one_ptm_mng_ptr->max_proteoform_mass_ == mng_ptr->max_proteoform_mass_ &&
         one_ptm_mng_ptr->filter_scale_ == mng_ptr->filter_scale_ &&
         one_ptm_mng_ptr->getIndexFilePara() == mng_ptr->getIndexFilePara();
}

bool linkOnePtmDiagIndex(const DiagFilterMngPtr& mng_ptr,
                         const OnePtmFilterMngPtr& one_ptm_mng_ptr,
                         const std::string& block_str) {
  namespace fs = std::filesystem;
  const std::vector<std::string>& one_names = one_ptm_mng_ptr->one_ptm_file_vec_;
  auto it = std::find(one_names.begin(), one_names.end(), "one_ptm_diag_index");
  if (it == one_names.end()) {
    return false;
  }
  std::string dir = indexDir(mng_ptr);
  fs::path src = dir + *it + one_ptm_mng_ptr->getIndexFilePara() + block_str;
  fs::path dst = dir + mng_ptr->multi_ptm_file_vec_[0] +
                 mng_ptr->getIndexFilePara() + block_str;
  std::error_code ec;
  if (!fs::is_regular_file(src, ec)) {
    return false;
  }
  if (fs::exists(dst, ec)) {
    // a rerun over an existing index: already the same file, or stale
    if (fs::equivalent(src, dst, ec)) {
      return true;
    }
    fs::remove(dst, ec);
  }
  fs::create_hard_link(src, dst, ec);
  if (!ec) {
    return true;
  }
  ec.clear();
  fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec);
  return !ec;
}

void geneDiagIndexFile(const ProteoformPtrVec& proteo_ptrs,
                       const DiagFilterMngPtr& mng_ptr, std::string block_str) {
  MassMatchPtr index_ptr = mass_match_factory::getPrmDiagMassMatchPtr(
      proteo_ptrs, mng_ptr->max_proteoform_mass_, mng_ptr->filter_scale_);

  std::string parameters = mng_ptr->getIndexFilePara();
  std::string dir_name = mng_ptr->prsm_para_ptr_->getOriDbName() + "_idx";
  std::string file_name =
      mng_ptr->multi_ptm_file_vec_[0] + parameters + block_str;
  {
    std::lock_guard<std::mutex> lock(serial_mutex);
    index_ptr->serializeMassMatch(file_name, dir_name);
  }
}

std::function<void()> geneIndexTask(int block_idx,
                                    const DiagFilterMngPtr& mng_ptr) {
  return [block_idx, mng_ptr]() {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string db_block_file_name =
        prsm_para_ptr->getSearchDbFileNameWithFolder() + "_" +
        std::to_string(block_idx);
    ProteoformPtrVec raw_forms =
        proteoform_factory::readFastaToProteoformPtrVec(
            db_block_file_name, prsm_para_ptr->getFixModPtrVec());
    std::string block_str = std::to_string(block_idx);
    geneDiagIndexFile(raw_forms, mng_ptr, block_str);
  };
}

void process(const DiagFilterMngPtr& mng_ptr,
             const OnePtmFilterMngPtr& one_ptm_mng_ptr) {
  std::string db_file_name =
      mng_ptr->prsm_para_ptr_->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  int block_num = db_block_ptr_vec.size();
  std::cout << "Generating multiple shift index files --- started" << std::endl;

  if (canReuseOnePtmDiagIndex(mng_ptr, one_ptm_mng_ptr)) {
    // the one-shift stage already wrote these files under another name
    for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
      std::cout << "Multiple shift index files - processing " << (i + 1)
                << " of " << block_num << " files." << std::endl;
      int block_idx = db_block_ptr_vec[i]->getBlockIdx();
      if (!linkOnePtmDiagIndex(mng_ptr, one_ptm_mng_ptr,
                               std::to_string(block_idx))) {
        LOG_WARN("Could not link the one-shift diagonal index of block "
                 << block_idx << "; building the multiple-shift index file.");
        geneIndexTask(block_idx, mng_ptr)();
      }
    }
    std::cout << "Generating multiple shift index files --- finished"
              << std::endl;
    return;
  }

  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(mng_ptr->thread_num_);
  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    std::cout << "Multiple shift index files - processing " << (i + 1) << " of "
              << block_num << " files." << std::endl;
    pool_ptr->enqueue(
        geneIndexTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr));
  }
  pool_ptr->shutDown();
  std::cout << "Generating multiple shift index files --- finished"
            << std::endl;
}

}  // namespace diag_index

}  // namespace toppic
