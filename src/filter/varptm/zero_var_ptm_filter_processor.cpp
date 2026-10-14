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
#include <memory>
#include <string>
#include <thread>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/logger.hpp"
#include "common/util/stage_profile.hpp"
#include "filter/massmatch/filter_spec.hpp"
#include "filter/oneptm/one_ptm_filter.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"
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

// The variable PTM and one-shift filters over the block, spectrum by
// spectrum. Both stages have one spectrum set per spectrum (precursor
// offset 0) in the same order unless the one-shift precursor offsets were
// changed; then the one-shift filter scores its spectra on its own. A
// spectrum whose mass lists are the same in both (the usual case) and fit
// below the one-shift column limit gives the one-shift filter the variable
// PTM filter's scores before the precursor bonus, so its four index
// scorings are saved; otherwise it scores the indexes with the windows cut
// at the limit, which gives its own files' scores.
void scoreVarOneBlock(const VarPtmFilterPtr& var_filter_ptr,
                      const OnePtmFilterPtr& one_filter_ptr, int block_idx,
                      const VarPtmFilterMngPtr& var_mng_ptr,
                      const OnePtmFilterMngPtr& one_mng_ptr,
                      const FilterSpecPtrVec& var_spectra,
                      const FilterSpecPtrVec& one_spectra,
                      SimplePrsmMemMerge* var_combiner,
                      SimplePrsmMemMerge* one_combiner) {
  PrsmParaPtr prsm_para_ptr = var_mng_ptr->prsm_para_ptr_;
  bool lockstep = one_spectra.size() == var_spectra.size();
  for (size_t i = 0; i < var_spectra.size(); i++) {
    const FilterSpecPtr& var_spec_ptr = var_spectra[i];
    bool one_done = false;
    std::function<void()> hook;
    if (lockstep) {
      const FilterSpecPtr& one_spec_ptr = one_spectra[i];
      hook = [&]() {
        if (one_filter_ptr->fitsColumns(one_spec_ptr) &&
            one_spec_ptr->getPrefMassErrors() ==
                var_spec_ptr->getPrefMassErrors() &&
            one_spec_ptr->getSuffMassErrors() ==
                var_spec_ptr->getSuffMassErrors()) {
          one_filter_ptr->computeBestMatch(
              one_spec_ptr, var_filter_ptr->getTermScores(),
              var_filter_ptr->getDiagScores(),
              var_filter_ptr->getRevTermScores(),
              var_filter_ptr->getRevDiagScores());
          one_done = true;
        }
      };
    }
    var_filter_ptr->computeBestMatch(var_spec_ptr, hook);
    var_ptm_filter_processor::addCands(var_filter_ptr, prsm_para_ptr,
                                       block_idx, static_cast<int>(i),
                                       var_combiner);
    var_ptm_filter_processor::reportProgress(var_mng_ptr, block_idx, i,
                                             var_spectra.size());
    if (lockstep) {
      if (!one_done) {
        one_filter_ptr->computeBestMatch(one_spectra[i]);
      }
      one_ptm_filter_processor::addCands(one_filter_ptr, prsm_para_ptr,
                                         block_idx, static_cast<int>(i),
                                         one_combiner);
      one_ptm_filter_processor::reportProgress(one_mng_ptr, block_idx, i,
                                               one_spectra.size());
    }
  }
  if (!lockstep) {
    one_ptm_filter_processor::scoreBlock(one_filter_ptr, block_idx,
                                         one_mng_ptr, one_spectra,
                                         one_combiner);
  }
}

std::function<void()> geneTask(int block_idx,
                               const ZeroPtmFilterMngPtr& zero_mng_ptr,
                               const VarPtmFilterMngPtr& var_mng_ptr,
                               const OnePtmFilterMngPtr& one_mng_ptr,
                               const FilterSpecPtrVec* zero_spectra,
                               const FilterSpecPtrVec* var_spectra,
                               const FilterSpecPtrVec* one_spectra,
                               SimplePrsmMemMerge* zero_combiner,
                               SimplePrsmMemMerge* var_combiner,
                               SimplePrsmMemMerge* one_combiner) {
  return [block_idx, zero_mng_ptr, var_mng_ptr, one_mng_ptr, zero_spectra,
          var_spectra, one_spectra, zero_combiner, var_combiner,
          one_combiner]() {
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
    OnePtmFilterPtr one_filter_ptr;
    if (one_mng_ptr != nullptr) {
      one_filter_ptr = std::make_shared<OnePtmFilter>(
          raw_forms, one_mng_ptr, zero_filter_ptr->getTermIndexPtr(),
          zero_filter_ptr->getDiagIndexPtr(),
          zero_filter_ptr->getRevTermIndexPtr(),
          zero_filter_ptr->getRevDiagIndexPtr());
    }
    zero_filter_ptr = nullptr;
    if (one_mng_ptr == nullptr) {
      var_ptm_filter_processor::scoreBlock(var_filter_ptr, block_idx,
                                           var_mng_ptr, *var_spectra,
                                           var_combiner);
    } else {
      scoreVarOneBlock(var_filter_ptr, one_filter_ptr, block_idx, var_mng_ptr,
                       one_mng_ptr, *var_spectra, *one_spectra, var_combiner,
                       one_combiner);
    }
  };
}

}  // namespace

void process(const ZeroPtmFilterMngPtr& zero_mng_ptr,
             const VarPtmFilterMngPtr& var_mng_ptr,
             const OnePtmFilterMngPtr& one_mng_ptr) {
  PrsmParaPtr prsm_para_ptr = zero_mng_ptr->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);
  int block_num = db_block_ptr_vec.size();
  int spec_num = msalign_util::getSpNum(sp_file_name);
  std::vector<int> sp_ids = msalign_util::readSpIds(sp_file_name);

  // the stages' candidate combiners, as their own processors build them
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

  std::unique_ptr<SimplePrsmMemMerge> one_combiner;
  if (one_mng_ptr != nullptr) {
    one_mng_ptr->n_spec_block_ = spec_num * block_num;
    one_mng_ptr->cnts_.resize(block_num, 0);
    std::vector<int> one_top_nums = {
        static_cast<int>(one_mng_ptr->comp_num_),
        static_cast<int>(one_mng_ptr->pref_suff_num_),
        static_cast<int>(one_mng_ptr->pref_suff_num_),
        static_cast<int>(one_mng_ptr->inte_num_)};
    one_combiner = std::make_unique<SimplePrsmMemMerge>(
        sp_file_name, one_mng_ptr->output_file_ext_, one_top_nums, sp_ids);
  }

  // the spectra of the stages, read and converted once; the three
  // readings (each 3-5 s on one thread) run at the same time
  FilterSpecPtrVec zero_spectra;
  FilterSpecPtrVec var_spectra;
  FilterSpecPtrVec one_spectra;
  {
    std::thread var_thread([&]() {
      var_spectra = var_ptm_filter_processor::readSpectra(var_mng_ptr);
    });
    std::thread one_thread([&]() {
      if (one_mng_ptr != nullptr) {
        one_spectra = one_ptm_filter_processor::readSpectra(one_mng_ptr);
      }
    });
    zero_spectra = zero_ptm_filter_processor::readSpectra(zero_mng_ptr);
    var_thread.join();
    one_thread.join();
  }

  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(zero_mng_ptr->thread_num_);
  LOG_DEBUG("thread num " << zero_mng_ptr->thread_num_);
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    pool_ptr->enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(),
                               zero_mng_ptr, var_mng_ptr, one_mng_ptr,
                               &zero_spectra, &var_spectra, &one_spectra,
                               &zero_combiner, &var_combiner,
                               one_combiner.get()));
  }
  pool_ptr->shutDown();
  stage_profile::report("zero_var_filter");
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
  if (one_combiner != nullptr) {
    std::cout << "One unexpected shift filtering - combining blocks started."
              << std::endl;
    one_combiner->write();
    std::cout << "One unexpected shift filtering - combining blocks finished."
              << std::endl;
  }
}

}  // namespace zero_var_ptm_filter_processor

}  // namespace toppic
