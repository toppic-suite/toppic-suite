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

#include "filter/varptm/var_ptm_filter.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <mutex>
#include <utility>
#include <vector>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/massmatch/index_load_gate.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "ms/factory/prm_ms_util.hpp"
#include "seq/prot_candidate.hpp"
#include "seq/proteoform_util.hpp"

namespace toppic {

// The first block of the stage is deserialized under the mutex (warming
// Boost.Serialization's lazy type-registration singletons on one thread);
// later blocks load in parallel — see the comment in zero_ptm_filter.cpp.
namespace {
std::mutex var_ptm_filter_mutex;
std::atomic<bool> var_ptm_filter_warmed{false};
}  // namespace

namespace {

// The precursor mass error shifted by each of shifts; shift 0 is skipped
// unless keep_zero (a zero shift is the zero-shift filter's, except for a
// proteoform with an N-terminal PTM, whose spectrum has the PTM mass taken
// off already).
std::vector<std::pair<int, int>> getShiftedMassErrors(
    std::pair<int, int> mass_error, const std::vector<int>& shifts,
    bool keep_zero) {
  std::vector<std::pair<int, int>> result;
  for (size_t i = 0; i < shifts.size(); i++) {
    if (shifts[i] == 0 && !keep_zero) {
      continue;
    }
    std::pair<int, int> m_e(mass_error.first - shifts[i], mass_error.second);
    result.push_back(m_e);
  }
  return result;
}

// The rows of an N-terminal side index whose first residue can carry a
// proteoform N-terminal PTM with the residues res_ptrs: the residue at the
// row's start is one of them and, in the N-terminal form index, the row's
// form is not an acetylated one (its N-terminus is blocked).
std::vector<char> compNTermRowMask(const MassMatchPtr& index_ptr,
                                   const std::vector<int>& row_starts,
                                   bool is_term_index,
                                   const ProteoformPtrVec& proteo_ptrs,
                                   const ProtModPtrVec& prot_mod_ptrs,
                                   const ResiduePtrVec& res_ptrs) {
  const std::vector<int>& row_proteo_ids = index_ptr->getRowProteoIds();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  std::vector<char> mask(row_proteo_ids.size(), 0);
  for (size_t r = 0; r < mask.size(); r++) {
    if (is_term_index) {
      bool acetylated = false;
      for (size_t m = 0; m < prot_mod_ptrs.size(); m++) {
        if (std::abs(prot_mod_ptrs[m]->getProtShift() - trunc_shifts[r]) <
            0.0001) {
          acetylated = prot_mod_ptrs[m]->isAcetylation();
          break;
        }
      }
      if (acetylated) {
        continue;
      }
    }
    const ResSeqPtr& seq_ptr =
        proteo_ptrs[row_proteo_ids[r]]->getResSeqPtr();
    int start = row_starts[r];
    if (start < 0 || start >= seq_ptr->getLen()) {
      continue;
    }
    const ResiduePtr& res_ptr = seq_ptr->getResiduePtr(start);
    for (size_t i = 0; i < res_ptrs.size(); i++) {
      if (res_ptr->isSame(res_ptrs[i])) {
        mask[r] = 1;
        break;
      }
    }
  }
  return mask;
}

// The candidates of one proteoform type from two rounds merged: a protein
// (by_range: a protein and a range) once, with its higher score, the first
// round's on a tie, and the best top_num in the order the finders report
// them: score descending, then protein id, start and end ascending.
ProtCandidatePtrVec mergeCands(const ProtCandidatePtrVec& first,
                               const ProtCandidatePtrVec& second,
                               size_t top_num, bool by_range) {
  ProtCandidatePtrVec all = first;
  for (size_t i = 0; i < second.size(); i++) {
    const ProtCandidatePtr& cand = second[i];
    bool found = false;
    for (size_t j = 0; j < all.size(); j++) {
      if (all[j]->getProteinId() == cand->getProteinId() &&
          (!by_range || (all[j]->getStartPos() == cand->getStartPos() &&
                         all[j]->getEndPos() == cand->getEndPos()))) {
        found = true;
        if (cand->getScore() > all[j]->getScore()) {
          all[j] = cand;
        }
        break;
      }
    }
    if (!found) {
      all.push_back(cand);
    }
  }
  std::stable_sort(all.begin(), all.end(),
                   [](const ProtCandidatePtr& a, const ProtCandidatePtr& b) {
                     if (a->getScore() != b->getScore()) {
                       return a->getScore() > b->getScore();
                     }
                     if (a->getProteinId() != b->getProteinId()) {
                       return a->getProteinId() < b->getProteinId();
                     }
                     if (a->getStartPos() != b->getStartPos()) {
                       return a->getStartPos() < b->getStartPos();
                     }
                     return a->getEndPos() < b->getEndPos();
                   });
  if (all.size() > top_num) {
    all.resize(top_num);
  }
  return all;
}

}  // namespace

VarPtmFilter::VarPtmFilter(const ProteoformPtrVec& proteo_ptrs,
                           const VarPtmFilterMngPtr& mng_ptr,
                           std::string block_str) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

  std::string index_dir = mng_ptr_->prsm_para_ptr_->getDbIndexDir();

  std::string parameters = mng_ptr->getIndexFilePara();
  std::string suffix = parameters + block_str;  // last part of file name

  // check if all index files for this ptm is present. if not, generate index
  // files again.
  bool index_files_exist = true;
  for (size_t t = 0; t < mng_ptr->var_ptm_file_vec_.size(); t++) {
    std::string file_name = mng_ptr->var_ptm_file_vec_[t] + suffix;
    if (!std::filesystem::exists(index_dir + file_util::getFileSeparator() +
                                 file_name)) {
      // if any of the index files for this ptm is missing
      index_files_exist = false;
      break;
    }
  }

  if (index_files_exist) {
    term_index_ptr_ = std::make_shared<MassMatch>();
    diag_index_ptr_ = std::make_shared<MassMatch>();
    rev_term_index_ptr_ = std::make_shared<MassMatch>();
    rev_diag_index_ptr_ = std::make_shared<MassMatch>();
    {
      LOG_DEBUG("Loading index files");
      std::unique_lock<std::mutex> lock(var_ptm_filter_mutex,
                                        std::defer_lock);
      if (!var_ptm_filter_warmed.load(std::memory_order_acquire)) {
        lock.lock();
        if (var_ptm_filter_warmed.load(std::memory_order_relaxed)) {
          lock.unlock();
        }
      }
      IndexLoadGate::Slot load_slot(index_load_gate);
      term_index_ptr_->deserializeMassMatch(
          mng_ptr->var_ptm_file_vec_[0] + suffix, index_dir);
      diag_index_ptr_->deserializeMassMatch(
          mng_ptr->var_ptm_file_vec_[1] + suffix, index_dir);
      rev_term_index_ptr_->deserializeMassMatch(
          mng_ptr->var_ptm_file_vec_[2] + suffix, index_dir);
      rev_diag_index_ptr_->deserializeMassMatch(
          mng_ptr->var_ptm_file_vec_[3] + suffix, index_dir);
      if (lock.owns_lock()) {
        var_ptm_filter_warmed.store(true, std::memory_order_release);
      }
    }
  } else {
    LOG_DEBUG("get shifts");
    std::vector<std::vector<double>> shift_2d =
        proteoform_util::getNTermShift2D(
            proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    std::vector<std::vector<double>> n_term_acet_2d =
        proteoform_util::getNTermAcet2D(
            proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    LOG_DEBUG("get shifts complete");
    // N-terminal indexes
    term_index_ptr_ = mass_match_factory::getPrmTermMassMatchPtr(
        proteo_ptrs, shift_2d, mng_ptr->max_proteoform_mass_,
        mng_ptr->filter_scale_);
    LOG_DEBUG("diag index");
    // Prm indexes
    diag_index_ptr_ = mass_match_factory::getPrmDiagMassMatchPtr(
        proteo_ptrs, mng_ptr->max_proteoform_mass_, mng_ptr->filter_scale_);
    std::vector<std::vector<double>> rev_shift_2d;
    std::vector<double> shift_1d(1, 0);
    for (size_t i = 0; i < proteo_ptrs.size(); i++) {
      rev_shift_2d.push_back(shift_1d);
    }
    // C-terminal indexes
    rev_term_index_ptr_ = mass_match_factory::getSrmTermMassMatchPtr(
        proteo_ptrs, rev_shift_2d, n_term_acet_2d,
        mng_ptr->max_proteoform_mass_, mng_ptr->filter_scale_);

    // To generate SRM indexes, n terminal acetylation shifts are added into the
    // SRM list.
    rev_diag_index_ptr_ = mass_match_factory::getSrmDiagMassMatchPtr(
        proteo_ptrs, n_term_acet_2d, mng_ptr->max_proteoform_mass_,
        mng_ptr->filter_scale_);
  }
  std::vector<int> proteo_lens;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    proteo_lens.push_back(proteo_ptrs[i]->getResSeqPtr()->getLen());
  }
  row_ranges_.init(proteo_lens, term_index_ptr_, diag_index_ptr_,
                   rev_term_index_ptr_, rev_diag_index_ptr_,
                   mng_ptr->prsm_para_ptr_->getProtModPtrVec());
  for (size_t k = 0; k < mng_ptr->n_term_shift_list_.size(); k++) {
    n_term_term_masks_.push_back(compNTermRowMask(
        term_index_ptr_, row_ranges_.term_starts_, true, proteo_ptrs,
        mng_ptr->prsm_para_ptr_->getProtModPtrVec(),
        mng_ptr->n_term_res_ptr_vec_2d_[k]));
    n_term_diag_masks_.push_back(compNTermRowMask(
        diag_index_ptr_, row_ranges_.diag_starts_, false, proteo_ptrs,
        mng_ptr->prsm_para_ptr_->getProtModPtrVec(),
        mng_ptr->n_term_res_ptr_vec_2d_[k]));
  }
}

FilterSpecPtr VarPtmFilter::geneFilterSpec(
    const PrmMsPtrVec& prm_ms_ptr_vec, const PrmMsPtrVec& srm_ms_ptr_vec,
    const VarPtmFilterMngPtr& mng_ptr) {
  PeakTolerancePtr tole_ptr =
      mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  std::vector<std::pair<int, int>> pref_mass_errors =
      prm_ms_util::getIntMassErrorList(prm_ms_ptr_vec, tole_ptr,
                                       mng_ptr->filter_scale_, true, false);
  std::vector<std::pair<int, int>> suff_mass_errors =
      prm_ms_util::getIntMassErrorList(srm_ms_ptr_vec, tole_ptr,
                                       mng_ptr->filter_scale_, false, true);
  std::pair<int, int> prec_minus_water_mass_error =
      prm_ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMassMinusWaterError(
          tole_ptr->getPpo(), mng_ptr->filter_scale_);
  return std::make_shared<FilterSpec>(
      prm_ms_ptr_vec[0]->getMsHeaderPtr(),
      static_cast<int>(prm_ms_ptr_vec.size()), std::move(pref_mass_errors),
      std::move(suff_mass_errors), prec_minus_water_mass_error);
}

void VarPtmFilter::scoreIndexes(
    const FilterSpecPtr& spec_ptr,
    const std::vector<std::pair<int, int>>& prec_minus_water_mass_errors) {
  const std::vector<std::pair<int, int>>& pref_mass_errors =
      spec_ptr->getPrefMassErrors();
  const std::vector<std::pair<int, int>>& suff_mass_errors =
      spec_ptr->getSuffMassErrors();

  term_scores_.assign(term_index_ptr_->getRowNum(), 0);
  term_prec_rows_.clear();
  term_index_ptr_->compMatchScores(pref_mass_errors,
                                   prec_minus_water_mass_errors, term_scores_,
                                   term_prec_rows_);

  diag_scores_.assign(diag_index_ptr_->getRowNum(), 0);
  diag_prec_rows_.clear();
  diag_index_ptr_->compMatchScores(pref_mass_errors,
                                   prec_minus_water_mass_errors, diag_scores_,
                                   diag_prec_rows_);

  rev_term_scores_.assign(rev_term_index_ptr_->getRowNum(), 0);
  rev_term_prec_rows_.clear();
  rev_term_index_ptr_->compMatchScores(suff_mass_errors,
                                       prec_minus_water_mass_errors,
                                       rev_term_scores_, rev_term_prec_rows_);

  rev_diag_scores_.assign(rev_diag_index_ptr_->getRowNum(), 0);
  rev_diag_prec_rows_.clear();
  rev_diag_index_ptr_->compMatchScores(suff_mass_errors,
                                       prec_minus_water_mass_errors,
                                       rev_diag_scores_, rev_diag_prec_rows_);

  // the bonus rows of each proteoform, collected once and shared by the four
  // proteoform types
  mass_match_util::compPrecMatchRows(term_index_ptr_, term_prec_rows_,
                                     term_rows_2d_);
  mass_match_util::compPrecMatchRows(diag_index_ptr_, diag_prec_rows_,
                                     diag_rows_2d_);
  mass_match_util::compPrecMatchRows(rev_term_index_ptr_, rev_term_prec_rows_,
                                     rev_term_rows_2d_);
  mass_match_util::compPrecMatchRows(rev_diag_index_ptr_, rev_diag_prec_rows_,
                                     rev_diag_rows_2d_);
}

std::array<ProtCandidatePtrVec, 4> VarPtmFilter::findTypeCands(
    const FilterSpecPtr& spec_ptr, std::vector<double>& ptm_shifts,
    const std::vector<char>* term_mask, const std::vector<char>* diag_mask) {
  PeakTolerancePtr tole_ptr =
      mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  double prec_minus_water_mass =
      spec_ptr->getHeaderPtr()->getFirstPrecMonoMassMinusWater();
  double prec_error_tole =
      spec_ptr->getHeaderPtr()->getFirstPrecErrorTolerance(tole_ptr->getPpo());

  // The candidates of a proteoform type (issue #348): with
  // range_candidates_, the best ranges, a protein once for each of its
  // ranges; otherwise the best proteins, each with the range of its best
  // pair.
  auto findCands = [&](const std::vector<short>& scores,
                       const std::vector<short>& rev_scores,
                       const mass_match_util::PrecMatchRows& rows_2d,
                       const mass_match_util::PrecMatchRows& rev_rows_2d,
                       const MassMatchPtr& index_ptr,
                       const MassMatchPtr& rev_index_ptr, bool term_start,
                       bool term_end, int threshold, int top_num) {
    mass_match_util::RangeRows range_rows =
        row_ranges_.getRangeRows(term_start, term_end, false, false);
    const std::vector<char>* start_row_mask =
        term_start ? term_mask : diag_mask;
    if (mng_ptr_->range_candidates_) {
      return mass_match_util::findVarPtmTopRanges(
          scores, rev_scores, rows_2d, rev_rows_2d, index_ptr, rev_index_ptr,
          range_rows, prec_minus_water_mass, prec_error_tole, ptm_shifts,
          threshold, top_num, start_row_mask);
    }
    return mass_match_util::findVarPtmTopProteins(
        scores, rev_scores, rows_2d, rev_rows_2d, index_ptr, rev_index_ptr,
        prec_minus_water_mass, prec_error_tole, ptm_shifts, threshold, top_num,
        &range_rows, start_row_mask);
  };

  std::array<ProtCandidatePtrVec, 4> cands;
  cands[0] = findCands(term_scores_, rev_term_scores_, term_rows_2d_,
                       rev_term_rows_2d_, term_index_ptr_, rev_term_index_ptr_,
                       true, true, mng_ptr_->comp_threshold_,
                       mng_ptr_->comp_num_);
  cands[1] = findCands(term_scores_, rev_diag_scores_, term_rows_2d_,
                       rev_diag_rows_2d_, term_index_ptr_, rev_diag_index_ptr_,
                       true, false, mng_ptr_->pref_suff_threshold_,
                       mng_ptr_->pref_suff_num_);
  cands[2] = findCands(diag_scores_, rev_term_scores_, diag_rows_2d_,
                       rev_term_rows_2d_, diag_index_ptr_, rev_term_index_ptr_,
                       false, true, mng_ptr_->pref_suff_threshold_,
                       mng_ptr_->pref_suff_num_);
  cands[3] = findCands(diag_scores_, rev_diag_scores_, diag_rows_2d_,
                       rev_diag_rows_2d_, diag_index_ptr_, rev_diag_index_ptr_,
                       false, false, mng_ptr_->internal_threshold_,
                       mng_ptr_->internal_num_);
  return cands;
}

void VarPtmFilter::computeBestMatch(const FilterSpecPtr& spec_ptr) {
  scoreIndexes(spec_ptr,
               getShiftedMassErrors(spec_ptr->getPrecMinusWaterMassError(),
                                    mng_ptr_->int_shift_list_, false));
  std::array<ProtCandidatePtrVec, 4> cands =
      findTypeCands(spec_ptr, mng_ptr_->shift_list_, nullptr, nullptr);

  // One more round for each proteoform N-terminal PTM mass m: the spectrum
  // with m taken off its N-terminal masses and its precursor mass, so that
  // the N-terminal rows match the fragments of a proteoform whose first
  // residue carries the PTM; the other variable PTMs are those of up to
  // var_ptm_num - 1, no other PTM included. Only the rows whose first
  // residue can carry the PTM start a pair.
  const FilterSpecPtrVec& n_term_specs = spec_ptr->getNTermModSpecs();
  const std::vector<size_t> top_nums = {
      mng_ptr_->comp_num_, mng_ptr_->pref_suff_num_, mng_ptr_->pref_suff_num_,
      mng_ptr_->internal_num_};
  for (size_t k = 0; k < n_term_specs.size(); k++) {
    const FilterSpecPtr& n_term_spec_ptr = n_term_specs[k];
    scoreIndexes(n_term_spec_ptr,
                 getShiftedMassErrors(
                     n_term_spec_ptr->getPrecMinusWaterMassError(),
                     mng_ptr_->n_term_other_int_shift_list_, true));
    std::array<ProtCandidatePtrVec, 4> round_cands = findTypeCands(
        n_term_spec_ptr, mng_ptr_->n_term_other_shift_list_,
        &n_term_term_masks_[k], &n_term_diag_masks_[k]);
    for (size_t t = 0; t < cands.size(); t++) {
      cands[t] = mergeCands(cands[t], round_cands[t], top_nums[t],
                            mng_ptr_->range_candidates_);
    }
  }

  int group_spec_num = spec_ptr->getGroupSpecNum();
  std::array<SimplePrsmPtrVec*, 4> match_ptrs = {
      &comp_match_ptrs_, &pref_match_ptrs_, &suff_match_ptrs_,
      &internal_match_ptrs_};
  for (size_t t = 0; t < cands.size(); t++) {
    match_ptrs[t]->clear();
    for (size_t i = 0; i < cands[t].size(); i++) {
      int id = cands[t][i]->getProteinId();
      match_ptrs[t]->push_back(std::make_shared<SimplePrsm>(
          spec_ptr->getHeaderPtr(), group_spec_num, proteo_ptrs_[id],
          cands[t][i]));
    }
  }
}

/*
void VarPtmFilter::computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec) {

  PeakTolerancePtr tole_ptr =
mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr(); bool pref =
true; std::vector<std::pair<int, int> > pref_mass_errors =
extend_ms_util::getExtendIntMassErrorList(ms_ptr_vec, pref,
mng_ptr_->filter_scale_); pref = false; std::vector<std::pair<int, int> >
suff_mass_errors = extend_ms_util::getExtendIntMassErrorList(ms_ptr_vec, pref,
mng_ptr_->filter_scale_); std::pair<int, int> prec_minus_water_mass_error =
ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMassMinusWaterError(tole_ptr->getPpo(),
                                                                        mng_ptr_->filter_scale_);
  std::vector<std::pair<int, int> > prec_minus_water_mass_errors
    = getShiftedMassErrors(prec_minus_water_mass_error,
mng_ptr_->int_shift_list_);

  int term_row_num = term_index_ptr_->getRowNum();
  std::vector<short> term_scores(term_row_num, 0);
  term_index_ptr_->compMatchScores(pref_mass_errors,
prec_minus_water_mass_errors, term_scores);

  int diag_row_num = diag_index_ptr_->getRowNum();
  std::vector<short> diag_scores(diag_row_num, 0);
  diag_index_ptr_->compMatchScores(pref_mass_errors,
prec_minus_water_mass_errors, diag_scores);

  int rev_term_row_num = rev_term_index_ptr_->getRowNum();
  std::vector<short> rev_term_scores(rev_term_row_num, 0);
  rev_term_index_ptr_->compMatchScores(suff_mass_errors,
prec_minus_water_mass_errors, rev_term_scores);

  int rev_diag_row_num = rev_diag_index_ptr_->getRowNum();
  std::vector<short> rev_diag_scores(rev_diag_row_num, 0);
  rev_diag_index_ptr_->compMatchScores(suff_mass_errors,
prec_minus_water_mass_errors, rev_diag_scores);

  int threshold = MassMatch::getPrecursorMatchScore() * 2 + 4;
  double prec_minus_water_mass =
ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMassMinusWater(); double
prec_error_tole =
ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecErrorTolerance(tole_ptr->getPpo());
  int group_spec_num = ms_ptr_vec.size();
  LOG_DEBUG("mass " << ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass());
  ProtCandidatePtrVec comp_prots
    = mass_match_util::findVarPtmTopProteins(term_scores, rev_term_scores,
                                             term_index_ptr_,
rev_term_index_ptr_, prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             threshold, mng_ptr_->comp_num_);
  comp_match_ptrs_.clear();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    int id = comp_prots[i]->getProteinId();
    comp_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id],
                                                            comp_prots[i]));
  }

  ProtCandidatePtrVec pref_prots
    = mass_match_util::findVarPtmTopProteins(term_scores, rev_diag_scores,
                                             term_index_ptr_,
rev_diag_index_ptr_, prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             threshold,
mng_ptr_->pref_suff_num_); pref_match_ptrs_.clear(); for (size_t i = 0; i <
pref_prots.size(); i++) { int id = pref_prots[i]->getProteinId();
    pref_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id],
                                                            pref_prots[i]));
  }

  ProtCandidatePtrVec suff_prots
    = mass_match_util::findVarPtmTopProteins(diag_scores, rev_term_scores,
                                             diag_index_ptr_,
rev_term_index_ptr_, prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             threshold,
mng_ptr_->pref_suff_num_); suff_match_ptrs_.clear(); for (size_t i = 0; i <
suff_prots.size(); i++) { int id = suff_prots[i]->getProteinId();
    suff_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id],
                                                            suff_prots[i]));
  }

  ProtCandidatePtrVec internal_prots
    = mass_match_util::findVarPtmTopProteins(diag_scores, rev_diag_scores,
                                             diag_index_ptr_,
rev_diag_index_ptr_, prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             threshold, mng_ptr_->inte_num_);
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    int id = internal_prots[i]->getProteinId();
    internal_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                                group_spec_num,
                                                                proteo_ptrs_[id],
                                                                internal_prots[i]));
  }
}
*/

} /* namespace toppic */
