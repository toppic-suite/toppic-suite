//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_TDGF_COMP_PVALUE_LOOKUP_TABLE_HPP_
#define TOPPIC_TDGF_COMP_PVALUE_LOOKUP_TABLE_HPP_

#include <vector>
#include <fstream>

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "tdgf/count_test_num.hpp"

namespace toppic {

class CompPValueLookupTable {
 public:
  explicit CompPValueLookupTable(TdgfMngPtr mng_ptr);

  bool inTable(const DeconvMsPtrVec &deconv_ms_ptr_vec, const PrsmPtrVec &prsm_ptrs);

  bool inTable(int peak_num, int match_frag_num, int unexpected_shift_num);

  void process(const DeconvMsPtrVec &deconv_ms_ptr_vec, PrsmPtrVec &prsm_ptrs, double ppo);

  double compProb(int peak_num, int match_frag_num, int unexpected_shift_num);

 private:
  void initTable();

  TdgfMngPtr mng_ptr_;

  CountTestNumPtr test_num_ptr_;

  std::ifstream input_;

  double ptm0_[48][20];

  double ptm1_[48][20];

  double ptm2_[48][20];
};

typedef std::shared_ptr<CompPValueLookupTable> CompPValueLookupTablePtr;

int getPeakIndex(int i);

int getFragIndex(int i);

// get x1, x2, y1, y2
std::vector<int> getFourIndex(int peak_num, int frag_num);

int getPeakNumFromIndex(int idx);

}  // namespace toppic

#endif
