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

#ifndef TOPPIC_SEARCH_GRAPHALIGN_GRAPH_POST_PROCESSOR_HPP_
#define TOPPIC_SEARCH_GRAPHALIGN_GRAPH_POST_PROCESSOR_HPP_

#include <string>
#include <unordered_map>
#include <vector>

#include "common/base/mod.hpp"
#include "common/base/prot_mod.hpp"
#include "search/graphalign/graph_align_mng.hpp"
#include "seq/mass_shift.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

// The protein N-terminal form of a PrSM, as GraphPostProcessor reads it from
// the alignment's start and first variable shift (chooseNTermForm).
struct GraphNTermForm {
  ProtModPtr prot_mod_ptr;
  // 1 when the proteoform starts one residue later than the alignment (an
  // NME_ACETYLATION read from an alignment that starts at the M), else 0
  int start_offset = 0;
  // whether the form's acetylation takes its mass from the first variable
  // shift; the mass left in that shift, and whether nothing is left
  bool takes_first_shift = false;
  double remain_mass = 0;
  bool drop_first_shift = false;
  // a truncation N-terminal modification (a proteoform N-terminal PTM of
  // the file, Position N-term) of a NONE proteoform starting at protein
  // residue 3 or later; nullptr when none
  ModPtr trunc_n_term_mod_ptr = nullptr;
};

class GraphPostProcessor {
 public:
  GraphPostProcessor(const GraphAlignMngPtr& mng_ptr,
                     const std::string& input_ext,
                     const std::string& output_ext)
      : mng_ptr_(mng_ptr), input_ext_(input_ext), output_ext_(output_ext) {}

  void process();

  void geneMassPtmMap();

 private:
  // the PTM combinations (mass_ptm_map_ entries) explaining the variable
  // shifts of one proteoform, at most max_known_mods_ PTMs in all
  std::vector<PtmPtrVec> choosePtmCombs(const MassShiftPtrVec& shift_vec);

  GraphAlignMngPtr mng_ptr_;

  std::string input_ext_;

  std::string output_ext_;

  std::unordered_map<int, PtmPtrVec> mass_ptm_map_;

  // the proteoform N-terminal PTMs of the file (Position N-term), one per
  // amino acid as mod_util::readModTxt lists them
  ModPtrVec n_term_mod_ptrs_;

  // the PTM combination of mass_ptm_map_ nearest to an integer mass, among
  // those of at most max_ptm_num PTMs; returns the integer mass error (-1
  // when there is none)
  int nearestPtmComb(int mass, int max_ptm_num, PtmPtrVec& ptm_vec);

  GraphNTermForm chooseNTermForm(const ProteoformPtr& form_ptr,
                                 const MassShiftPtr& first_shift_ptr);
};

using GraphPostProcessorPtr = std::shared_ptr<GraphPostProcessor>;

}  // namespace toppic

#endif
