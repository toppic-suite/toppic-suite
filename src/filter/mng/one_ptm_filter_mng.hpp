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

#ifndef TOPPIC_FILTER_MNG_ONE_PTM_FILTER_MNG_HPP_
#define TOPPIC_FILTER_MNG_ONE_PTM_FILTER_MNG_HPP_

#include <mutex>
#include <set>

#include "filter/mng/ptm_diag_index.hpp"
#include "para/prsm_para.hpp"

namespace toppic {

class OnePtmFilterMng {
 public:
  // for index generation
  OnePtmFilterMng(const PrsmParaPtr& prsm_para_ptr,
                  const std::string& index_file_para,
                  const std::string& output_file_ext, int thread_num);

  // for one shift filter in TopPIC
  OnePtmFilterMng(const PrsmParaPtr& prsm_para_ptr,
                  const std::string& index_file_para,
                  const std::string& output_file_ext, int thread_num,
                  double min_shift, double max_shift);

  // for One shift filter in TopMG
  OnePtmFilterMng(const PrsmParaPtr& prsm_para_ptr,
                  const std::string& index_file_para,
                  const std::string& output_file_ext, int thread_num,
                  double min_shift, double max_shift,
                  const std::string& residueModFileName, int var_num);

  std::string getIndexFilePara() { return index_file_para_; }

  PrsmParaPtr prsm_para_ptr_;

  // parameters for fast filteration; the column limits of the whole
  // one-shift index family, which the multiple-shift filter shares through
  // its diagonal file (see ptm_diag_index.hpp)
  int max_proteoform_mass_ = ptm_diag_index::kMaxProteoformMass;

  // Candidate protein number for each spectrum
  unsigned int comp_num_ = 10;
  unsigned int pref_suff_num_ = 10;
  unsigned int inte_num_ = 10;

  // reported n term or c term shift number
  unsigned int comp_cand_shift_num_ = 1;
  unsigned int cand_shift_num_ = 10;
  int filter_scale_ = ptm_diag_index::kFilterScale;

  int threshold_ = 4;

  std::string index_file_para_;

  std::string output_file_ext_;

  int thread_num_;

  double min_shift_ = -500;
  double max_shift_ = 500;

  std::string residueModFileName_ = "";

  int var_num_ = 0;

  // true: a candidate is a protein and a sequence range, so a protein can
  // take several of the candidate numbers above with different ranges (issue
  // #348); false (the default): each protein keeps its best range and the
  // shifts of its best pairs. Here the pairs of the shift window give a
  // protein many ranges that share one strong side, which took 29% to 55% of
  // the slots on the 480 data set and made the one-shift search worse on
  // five times as many spectra as it made it better (toppic_test experiment
  // 86).
  bool range_candidates_ = false;

  // The spectra the variable-PTM variants are filtered for (var_num_ > 0),
  // by spectrum id. A topmg run filters every spectrum and leaves
  // search_all_ true; toppic's --topmg filters only the spectra its graph
  // alignment searches: with pipeline the ones its earlier stages have not
  // identified, with topmg-only all of them.
  bool search_all_ = true;

  std::set<int> search_spec_ids_;

  int n_spec_block_ = 0;

  std::mutex mutex_;

  std::vector<int> cnts_;

  // file name vector: N-terminal, PRM diagonal (the file the multiple-shift
  // filter reads as well), C-terminal, SRM diagonal
  std::vector<std::string> one_ptm_file_vec_{
      "one_ptm_term_index", ptm_diag_index::kFileName,
      "one_ptm_rev_term_index", "one_ptm_rev_diag_index"};
};

using OnePtmFilterMngPtr = std::shared_ptr<OnePtmFilterMng>;

}  // namespace toppic

#endif /* ONE_PTM_FILTER_MNG_HPP_ */
