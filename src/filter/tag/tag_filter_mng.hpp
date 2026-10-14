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

#ifndef TOPPIC_FILTER_TAG_TAG_FILTER_MNG_HPP_
#define TOPPIC_FILTER_TAG_TAG_FILTER_MNG_HPP_

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "para/prsm_para.hpp"

namespace toppic {

// The settings of the tag filter (filter/tag, toppic's hidden --tag-filter)
// for one identification stage. The stage decides the precursor error list
// of the spectrum sets, what a candidate carries (shifts, a sequence range
// or the protein alone) and the files written (one per proteoform type or
// one in all).
class TagFilterMng {
 public:
  enum class Stage { kVar, kOne, kMulti, kTopmg };

  TagFilterMng(const PrsmParaPtr& prsm_para_ptr, Stage stage, int thread_num,
               const std::string& output_file_ext);

  // the option names of the stages, in the order of Stage
  static const std::vector<std::string>& getStageNames();

  static const std::string& getStageName(Stage stage);

  // The stages named by a --tag-filter value: a comma-separated list of
  // stage names, or "all". Returns false when a name is not a stage.
  static bool parseStages(const std::string& value, std::set<Stage>* stages);

  // the precursor error list of the stage's spectrum sets
  const std::vector<double>& getPrecErrorVec() const;

  // whether the candidates of a stage are written per proteoform type
  bool writeTypeFiles() const {
    return stage_ == Stage::kVar || stage_ == Stage::kOne;
  }

  // whether a spectrum (by id) is scored
  bool searchSpectrum(int spec_id) const;

  PrsmParaPtr prsm_para_ptr_;
  Stage stage_;
  int thread_num_ = 1;
  std::string output_file_ext_;
  // --tag-filter-length: the residues of a tag (at least TagIndex::K)
  int tag_len_ = 4;
  // --tag-filter-cand-num: the proteins kept per spectrum set
  int cand_num_ = 10;
  // The least score of a candidate, and the least support of a second
  // diagonal that counts towards the score: tag_len_ + 3 peaks, a ladder of
  // tag_len_ + 2 residues (three overlapping tags) or two tags apart on one
  // diagonal (0 = that default). A protein with one tag, or a ladder one
  // residue longer, is not a candidate: a tag of four residues occurs in
  // hundreds of proteins by chance and a ladder of five in several, so
  // those tiers hold the true protein among many chance ones.
  int min_score_ = 0;

  int getMinScore() const { return min_score_ > 0 ? min_score_ : tag_len_ + 3; }
  // the distinct diagonal shifts a one-shift candidate carries at most
  int max_shift_num_ = 10;
  // two diagonal shifts closer than this are one shift
  double shift_merge_tole_ = 0.01;
  // var stage: the variable PTM file and number (-b, --num-variable-ptms)
  std::string var_ptm_file_name_;
  int var_ptm_num_ = 0;
  // The spectra scored: every spectrum when search_all_, else those of
  // search_spec_ids_ (the hosted TopMG stage), without those of
  // skip_spec_ids_ (the multiple-shift stage) in both cases.
  bool search_all_ = true;
  std::set<int> search_spec_ids_;
  std::set<int> skip_spec_ids_;
};

using TagFilterMngPtr = std::shared_ptr<TagFilterMng>;

}  // namespace toppic

#endif
