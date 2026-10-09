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

#ifndef TOPPIC_STAT_TDGF_EVALUE_PROCESSOR_HPP_
#define TOPPIC_STAT_TDGF_EVALUE_PROCESSOR_HPP_

#include "prsm/prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "stat/tdgf/comp_pvalue_array.hpp"
#include "stat/tdgf/tdgf_mng.hpp"

namespace toppic {

class EValueProcessor {
 public:
  EValueProcessor(const TdgfMngPtr& mng_ptr) : mng_ptr_(mng_ptr) {};

  void init();

  // Builds the CountTestNum of mng_ptr's parameters into the cache that
  // init() reads (a database scan of about 9 s on one thread), so that the
  // scan can run in a background thread during the filter stages instead
  // of at the start of the E-value stage (toppic_test experiment 150).
  static void prepareCountTestNum(const TdgfMngPtr& mng_ptr);

  void process(bool is_separate);

  void compEvalues(const SpectrumSetPtr& spectrum_set_ptr,
                   PrsmPtrVec& sele_prsm_ptrs, bool is_separate);

  // compEvalues through a per-worker-thread CompPValueArray (the same
  // pattern the process() thread pool uses), for callers that run one
  // spectrum's E-values per pool task. Safe to call concurrently for
  // disjoint PrSM sets.
  void compEvaluesThreadLocal(const SpectrumSetPtr& spectrum_set_ptr,
                              PrsmPtrVec& sele_prsm_ptrs, bool is_separate);

 private:
  TdgfMngPtr mng_ptr_;

  CountTestNumPtr test_num_ptr_;

  CompPValueArrayPtr comp_pvalue_ptr_;

  bool checkPrsms(const PrsmPtrVec& prsm_ptrs);

  void fixupEvalues(PrsmPtrVec& sele_prsm_ptrs);
};

using EValueProcessorPtr = std::shared_ptr<EValueProcessor>;

}  // namespace toppic

#endif
