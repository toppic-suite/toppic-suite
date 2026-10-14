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

#ifndef TOPPIC_FILTER_TAG_TAG_FILTER_HPP_
#define TOPPIC_FILTER_TAG_TAG_FILTER_HPP_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "filter/tag/tag_filter_mng.hpp"
#include "filter/tag/tag_index.hpp"
#include "ms/spec/prm_peak.hpp"
#include "ms/spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

// The tag filter of one spectrum set (toppic's hidden --tag-filter).
//
// Tags: the PRM peaks of each base type (ORIGINAL, the N-terminal PRM of a
// fragment mass; REVERSED, the precursor minus the mass), peaks within
// their tolerances of each other merged (one fragment at several charge
// states is listed once per charge in the msalign file). Two peaks whose
// difference is a residue mass within the sum of their strict tolerances
// are joined by an edge labelled with the residue (every residue that
// fits). A tag is a path of tag_len_ edges, read in increasing mass, so it
// reads the protein from the N-terminus whichever ion series it comes
// from; its string has I written as L.
//
// Candidates: every occurrence of every tag in the database (TagIndex). An
// occurrence lies on the diagonal s = (PRM mass of the tag's first peak) -
// (the protein's prefix mass before the tag), and the occurrences of a
// protein whose diagonals agree within their peaks' tolerances form one
// diagonal. A protein's score is the number of distinct PRM peaks covered
// by the tags of its best diagonal, plus those of its second diagonal when
// two or more tags lie on it (the two parts of a proteoform with one
// modification), so that the chance hits spread over the diagonals of a
// long protein do not add up. The cand_num_ proteins of the highest score
// (ties: the shorter protein, then the name) with a score of at least
// TagFilterMng::getMinScore (a ladder of tag_len_ + 2 residues or two tags
// apart on one diagonal) are the candidates. What their simple PrSMs
// carry depends on the stage (TagFilterMng::Stage):
//   one    the diagonals' shifts as n_term_shifts and, as c = P - protein
//          mass - s, as c_term_shifts, so that the one-shift search aligns
//          on the tags' diagonals; the same candidate for every proteoform
//          type;
//   var    the proteoforms containing a tag whose mass matches the
//          precursor up to variable PTMs: a sequence range, the protein
//          N-terminal form's shift (complete and prefix) or a truncation
//          start (suffix and internal), a truncation N-terminal
//          modification when one is needed, by type;
//   multi, topmg   the protein alone.
class TagFilter {
 public:
  TagFilter(const TagIndexPtr& index_ptr, const TagFilterMngPtr& mng_ptr);

  // statistics of one spectrum set for the processor's summary
  struct Stats {
    size_t tag_num = 0;
    size_t occurrence_num = 0;
    size_t prot_num = 0;
  };

  // The candidates of a spectrum set, by proteoform type (COMPLETE, PREFIX,
  // SUFFIX, INTERNAL); the multi and topmg stages put theirs in [0].
  std::vector<SimplePrsmPtrVec> process(
      const SpectrumSetPtr& spec_set_ptr) const;

  std::vector<SimplePrsmPtrVec> process(const SpectrumSetPtr& spec_set_ptr,
                                        Stats* stats) const;

 private:
  struct TagResidue {
    char letter;
    double mass;
  };

  // a tag: its string and the ids of its tag_len_ + 1 peaks, in increasing
  // mass
  struct Tag {
    std::string str;
    std::vector<int> peak_ids;
  };

  // an occurrence of a tag in a protein
  struct Occurrence {
    // the 0-based position of the tag's first residue
    int pos;
    // the PRM mass of the tag's first peak and its N-terminal tolerance
    double first_mass;
    double first_tole;
    // the diagonal shift first_mass - prefix mass before pos
    double shift;
    std::vector<int> peak_ids;
  };

  // a diagonal of a protein: the mean shift of its occurrences and the
  // number of distinct peaks they cover
  struct Diagonal {
    double shift;
    int support;
  };

  struct ProtHit {
    std::vector<Occurrence> occs;
    std::vector<Diagonal> diagonals;
    int score = 0;
  };

  struct Candidate {
    int prot_idx;
    int score;
    std::vector<Occurrence> occs;
    // the diagonals, the best supported first
    std::vector<Diagonal> diagonals;
  };

  void initResidues();

  void initVarPtmShifts();

  // the tags of the peaks (sorted by position, ids set)
  void geneTags(const PrmPeakPtrVec& peaks, std::vector<Tag>* tags) const;

  void extendTag(const std::vector<std::vector<std::pair<int, char>>>& edges,
                 std::string* str, std::vector<int>* peak_ids,
                 std::vector<Tag>* tags) const;

  // the diagonals of a protein's occurrences (sorted by shift afterwards)
  // and its score
  void scoreHit(ProtHit* hit) const;

  // the candidates of the peaks: tags, their occurrences, the top proteins
  std::vector<Candidate> findCandidates(const PrmPeakPtrVec& peaks,
                                        Stats* stats) const;

  void geneOneShiftPrsms(const Candidate& cand, const MsHeaderPtr& header_ptr,
                         int group_spec_num,
                         std::vector<SimplePrsmPtrVec>* results) const;

  void geneVarPtmPrsms(const Candidate& cand, const MsHeaderPtr& header_ptr,
                       int group_spec_num,
                       std::vector<SimplePrsmPtrVec>* results) const;

  void geneProtPrsm(const Candidate& cand, const MsHeaderPtr& header_ptr,
                    int group_spec_num,
                    std::vector<SimplePrsmPtrVec>* results) const;

  TagIndexPtr index_ptr_;
  TagFilterMngPtr mng_ptr_;
  // the residues of the tags with the fixed modifications, by mass
  std::vector<TagResidue> residues_;
  double min_res_mass_ = 0;
  double max_res_mass_ = 0;
  // var stage: the sums of 0 to var_ptm_num_ variable PTM masses, the sums
  // of one PTM fewer (with a truncation N-terminal modification), and the
  // truncation N-terminal modification masses (Position N-term)
  std::vector<double> var_shift_sums_;
  std::vector<double> var_shift_sums_fewer_;
  std::vector<double> n_term_mod_shifts_;
};

using TagFilterPtr = std::shared_ptr<TagFilter>;

}  // namespace toppic

#endif
