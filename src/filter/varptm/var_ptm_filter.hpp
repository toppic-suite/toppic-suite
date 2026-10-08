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

#ifndef TOPPIC_FILTER_VARPTM_VAR_PTM_FILTER_HPP_
#define TOPPIC_FILTER_VARPTM_VAR_PTM_FILTER_HPP_

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "filter/massmatch/filter_spec.hpp"
#include "filter/massmatch/mass_match.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"
#include "ms/spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"

namespace toppic {

class VarPtmFilter {
 public:
  VarPtmFilter(const ProteoformPtrVec& proteo_ptrs,
               const VarPtmFilterMngPtr& mng_ptr, std::string block_str);

  // The filter on indexes another filter has loaded: the zero-shift
  // filter's four indexes of the same block are this filter's (the same
  // files), so the combined zero-shift / variable PTM stage
  // (zero_var_ptm_filter_processor) loads each block once (toppic_test
  // experiment 150).
  VarPtmFilter(const ProteoformPtrVec& proteo_ptrs,
               const VarPtmFilterMngPtr& mng_ptr, MassMatchPtr term_index_ptr,
               MassMatchPtr diag_index_ptr, MassMatchPtr rev_term_index_ptr,
               MassMatchPtr rev_diag_index_ptr);

  // the block independent part of a spectrum set, computed once per stage
  static FilterSpecPtr geneFilterSpec(const PrmMsPtrVec& prm_ms_ptr_vec,
                                     const PrmMsPtrVec& srm_ms_ptr_vec,
                                     const VarPtmFilterMngPtr& mng_ptr);

  // The candidates of a spectrum set. With truncation N-terminal
  // modifications (variable PTMs of the proteoform N-terminus), the
  // spectrum is scored once more for each of their masses (the spectra of
  // spec_ptr->getNTermModSpecs()), and each of these rounds has candidates
  // of its own (getNTermMatchPtrs, carrying the round's mass in
  // SimplePrsm::getTruncNTermModShifts), so that the many weak pairs a round
  // with the precursor mass less the PTM mass finds do not displace the
  // candidates of the first round.
  void computeBestMatch(const FilterSpecPtr& spec_ptr);

  SimplePrsmPtrVec getCompMatchPtrs() { return comp_match_ptrs_; }
  SimplePrsmPtrVec getPrefMatchPtrs() { return pref_match_ptrs_; }
  SimplePrsmPtrVec getSuffMatchPtrs() { return suff_match_ptrs_; }
  SimplePrsmPtrVec getInternalMatchPtrs() { return internal_match_ptrs_; }

  // the candidates of the N-terminal round k of the last spectrum set, of
  // proteoform type type_idx (0 complete, 1 prefix, 2 suffix, 3 internal)
  const SimplePrsmPtrVec& getNTermMatchPtrs(size_t k, int type_idx) {
    return n_term_match_ptrs_[k][type_idx];
  }

  size_t getNTermRoundNum() { return n_term_match_ptrs_.size(); }

  MassMatchPtr getDiagIndexPtr() { return diag_index_ptr_; }
  MassMatchPtr getRevDiagIndexPtr() { return rev_diag_index_ptr_; }
  MassMatchPtr getTermIndexPtr() { return term_index_ptr_; }
  MassMatchPtr getRevTermIndexPtr() { return rev_term_index_ptr_; }

 private:
  VarPtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;

  // Score arrays of the four indexes, the rows of each that received the
  // precursor bonus, and those rows grouped by proteoform, which the four
  // proteoform types share. They are members so that they are allocated once
  // per database block instead of once per spectrum.
  std::vector<short> term_scores_;
  std::vector<short> diag_scores_;
  std::vector<short> rev_term_scores_;
  std::vector<short> rev_diag_scores_;
  std::vector<int> term_prec_rows_;
  std::vector<int> diag_prec_rows_;
  std::vector<int> rev_term_prec_rows_;
  std::vector<int> rev_diag_prec_rows_;
  mass_match_util::PrecMatchRows term_rows_2d_;
  mass_match_util::PrecMatchRows diag_rows_2d_;
  mass_match_util::PrecMatchRows rev_term_rows_2d_;
  mass_match_util::PrecMatchRows rev_diag_rows_2d_;

  MassMatchPtr diag_index_ptr_;
  MassMatchPtr rev_diag_index_ptr_;
  MassMatchPtr term_index_ptr_;
  MassMatchPtr rev_term_index_ptr_;

  // the sequence ranges of the index rows
  mass_match_util::IndexRowRanges row_ranges_;

  // the row ranges and the N-terminal masks, from the four indexes
  void initRows();

  // For each truncation N-terminal modification mass of the manager, the
  // rows of the diagonal index whose first residue can carry the PTM: 1 if
  // the row starts at the third residue of the protein or later and the
  // amino acid there is one of the PTM's.
  std::vector<std::vector<char>> n_term_diag_masks_;

  // Score the spectrum against the four indexes (without score_term, all
  // but the N-terminal form index), with the precursor bonus for the
  // precursor masses of prec_minus_water_mass_errors.
  void scoreIndexes(const FilterSpecPtr& spec_ptr,
                    const std::vector<std::pair<int, int>>&
                        prec_minus_water_mass_errors,
                    bool score_term);

  // The candidates of the four proteoform types (complete, prefix, suffix,
  // internal) from the scores of scoreIndexes: the pairs whose range mass
  // plus one of ptm_shifts is the precursor mass of spec_ptr. With
  // diag_mask (a truncation N-terminal modification round), only the
  // suffix and internal candidates, whose diagonal index rows are
  // restricted to diag_mask.
  std::array<ProtCandidatePtrVec, 4> findTypeCands(
      const FilterSpecPtr& spec_ptr, std::vector<double>& ptm_shifts,
      const std::vector<char>* diag_mask);

  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
  std::vector<std::array<SimplePrsmPtrVec, 4>> n_term_match_ptrs_;

  // which db_block it is. Used when merging the output files by order.
  std::string block_str;
};

using VarPtmFilterPtr = std::shared_ptr<VarPtmFilter>;

} /* namespace toppic */

#endif
