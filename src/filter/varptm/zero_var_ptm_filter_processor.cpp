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

#include "filter/varptm/zero_var_ptm_filter_processor.hpp"

#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <thread>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/logger.hpp"
#include "filter/massmatch/filter_spec.hpp"
#include "filter/varptm/var_ptm_filter.hpp"
#include "filter/varptm/var_ptm_filter_processor.hpp"
#include "filter/zeroptm/zero_ptm_filter.hpp"
#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/simple_prsm_mem_merge.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace zero_var_ptm_filter_processor {

namespace {

std::function<void()> geneTask(int block_idx,
                               const ZeroPtmFilterMngPtr& zero_mng_ptr,
                               const VarPtmFilterMngPtr& var_mng_ptr,
                               const FilterSpecPtrVec* zero_spectra,
                               const FilterSpecPtrVec* var_spectra,
                               SimplePrsmMemMerge* zero_combiner,
                               SimplePrsmMemMerge* var_combiner) {
  return [block_idx, zero_mng_ptr, var_mng_ptr, zero_spectra, var_spectra,
          zero_combiner, var_combiner]() {
    PrsmParaPtr prsm_para_ptr = zero_mng_ptr->prsm_para_ptr_;
    std::string db_block_file_name =
        prsm_para_ptr->getSearchDbFileNameWithFolder() + "_" +
        std::to_string(block_idx);
    ProteoformPtrVec raw_forms =
        proteoform_factory::readFastaToProteoformPtrVec(
            db_block_file_name, prsm_para_ptr->getFixModPtrVec());
    // the zero-shift filter loads the block's four indexes
    ZeroPtmFilterPtr zero_filter_ptr = std::make_shared<ZeroPtmFilter>(
        raw_forms, zero_mng_ptr, std::to_string(block_idx));
    zero_ptm_filter_processor::scoreBlock(zero_filter_ptr, block_idx,
                                          zero_mng_ptr, *zero_spectra,
                                          zero_combiner);
    // the variable PTM filter scores on the same indexes
    VarPtmFilterPtr var_filter_ptr = std::make_shared<VarPtmFilter>(
        raw_forms, var_mng_ptr, zero_filter_ptr->getTermIndexPtr(),
        zero_filter_ptr->getDiagIndexPtr(),
        zero_filter_ptr->getRevTermIndexPtr(),
        zero_filter_ptr->getRevDiagIndexPtr());
    zero_filter_ptr = nullptr;
    var_ptm_filter_processor::scoreBlock(var_filter_ptr, block_idx,
                                         var_mng_ptr, *var_spectra,
                                         var_combiner);
  };
}

}  // namespace

void process(const ZeroPtmFilterMngPtr& zero_mng_ptr,
             const VarPtmFilterMngPtr& var_mng_ptr) {
  PrsmParaPtr prsm_para_ptr = zero_mng_ptr->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);
  int block_num = db_block_ptr_vec.size();
  int spec_num = msalign_util::getSpNum(sp_file_name);
  std::vector<int> sp_ids = msalign_util::readSpIds(sp_file_name);

  // the two stages' candidate combiners, as their own processors build them
  zero_mng_ptr->n_spec_block_ = spec_num * block_num;
  zero_mng_ptr->cnts_.resize(block_num, 0);
  std::vector<int> zero_top_nums = {
      static_cast<int>(zero_mng_ptr->comp_num_),
      static_cast<int>(zero_mng_ptr->pref_suff_num_),
      static_cast<int>(zero_mng_ptr->pref_suff_num_),
      static_cast<int>(zero_mng_ptr->inte_num_)};
  SimplePrsmMemMerge zero_combiner(sp_file_name, zero_mng_ptr->output_file_ext_,
                                   zero_top_nums, sp_ids);

  var_mng_ptr->n_spec_block_ = spec_num * block_num;
  var_mng_ptr->cnts_.resize(block_num, 0);
  std::vector<int> var_top_nums = {
      static_cast<int>(var_mng_ptr->comp_num_),
      static_cast<int>(var_mng_ptr->pref_suff_num_),
      static_cast<int>(var_mng_ptr->pref_suff_num_),
      static_cast<int>(var_mng_ptr->internal_num_)};
  SimplePrsmMemMerge var_combiner(
      sp_file_name, var_mng_ptr->output_file_ext_, var_top_nums, sp_ids,
      1 + static_cast<int>(var_mng_ptr->n_term_shift_list_.size()));

  // the spectra of both stages, read and converted once
  FilterSpecPtrVec zero_spectra =
      zero_ptm_filter_processor::readSpectra(zero_mng_ptr);
  FilterSpecPtrVec var_spectra =
      var_ptm_filter_processor::readSpectra(var_mng_ptr);

  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(zero_mng_ptr->thread_num_);
  LOG_DEBUG("thread num " << zero_mng_ptr->thread_num_);
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    pool_ptr->enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(),
                               zero_mng_ptr, var_mng_ptr, &zero_spectra,
                               &var_spectra, &zero_combiner, &var_combiner));
  }
  pool_ptr->shutDown();
  std::cout << std::endl;
  std::cout << "Zero unexpected shift filtering - combining blocks started."
            << std::endl;
  zero_combiner.write();
  std::cout << "Zero unexpected shift filtering - combining blocks finished."
            << std::endl;
  std::cout << "Variable PTM filtering - combining blocks started."
            << std::endl;
  var_combiner.write();
  std::cout << "Variable PTM filtering - combining blocks finished."
            << std::endl;
}

}  // namespace zero_var_ptm_filter_processor

}  // namespace toppic
