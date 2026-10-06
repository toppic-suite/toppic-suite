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

#include "filter/oneptm/one_ptm_filter.hpp"

#include <atomic>
#include <filesystem>
#include <iostream>
#include <mutex>

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
// On the 480 benchmark this stage read 84 s of index one block at a time
// through the old always-held lock (toppic_test experiment 66).
namespace {
std::mutex one_ptm_filter_mutex;
std::atomic<bool> one_ptm_filter_warmed{false};
}  // namespace
OnePtmFilter::OnePtmFilter(const ProteoformPtrVec& proteo_ptrs,
                           const OnePtmFilterMngPtr& mng_ptr,
                           std::string block_str) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

  std::string index_dir = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx";

  std::string parameters = mng_ptr->getIndexFilePara();
  std::string suffix = parameters + block_str;

  bool index_files_exist = true;

  for (size_t t = 0; t < mng_ptr->one_ptm_file_vec_.size(); t++) {
    std::string file_name = mng_ptr->one_ptm_file_vec_[t] + suffix;
    if (!std::filesystem::exists(index_dir + file_util::getFileSeparator() +
                                 file_name)) {
      // if any of the index files for this ptm is missing
      index_files_exist = false;
      // the diagonal file's earlier names
      mass_match_util::warnOldIndexLayout(index_dir, suffix);
      break;
    }
  }

  if (index_files_exist) {  // if exists
    term_index_ptr_ = std::make_shared<MassMatch>();
    diag_index_ptr_ = std::make_shared<MassMatch>();
    rev_term_index_ptr_ = std::make_shared<MassMatch>();
    rev_diag_index_ptr_ = std::make_shared<MassMatch>();

    {
      std::unique_lock<std::mutex> lock(one_ptm_filter_mutex,
                                        std::defer_lock);
      if (!one_ptm_filter_warmed.load(std::memory_order_acquire)) {
        lock.lock();
        if (one_ptm_filter_warmed.load(std::memory_order_relaxed)) {
          lock.unlock();
        }
      }
      IndexLoadGate::Slot load_slot(index_load_gate);
      // complete file name is a ptm type from topindexfile vector + parameters
      // + db block id (block_str)
      term_index_ptr_->deserializeMassMatch(
          mng_ptr->one_ptm_file_vec_[0] + suffix, index_dir);
      diag_index_ptr_->deserializeMassMatch(
          mng_ptr->one_ptm_file_vec_[1] + suffix, index_dir);
      rev_term_index_ptr_->deserializeMassMatch(
          mng_ptr->one_ptm_file_vec_[2] + suffix, index_dir);
      rev_diag_index_ptr_->deserializeMassMatch(
          mng_ptr->one_ptm_file_vec_[3] + suffix, index_dir);
      if (lock.owns_lock()) {
        one_ptm_filter_warmed.store(true, std::memory_order_release);
      }
    }
    // the files must have been built at this filter's column limits (the
    // diagonal one is shared with the multiple-shift filter)
    MassMatchPtr read_ptrs[] = {term_index_ptr_, diag_index_ptr_,
                                rev_term_index_ptr_, rev_diag_index_ptr_};
    for (size_t t = 0; t < 4; t++) {
      if (!mass_match_util::hasFilterLimits(read_ptrs[t],
                                            mng_ptr->max_proteoform_mass_,
                                            mng_ptr->filter_scale_)) {
        mass_match_util::warnIndexLimits(
            mng_ptr->one_ptm_file_vec_[t] + suffix, read_ptrs[t],
            mng_ptr->max_proteoform_mass_, mng_ptr->filter_scale_);
        index_files_exist = false;
        break;
      }
    }
    if (!index_files_exist) {
      term_index_ptr_ = nullptr;
      diag_index_ptr_ = nullptr;
      rev_term_index_ptr_ = nullptr;
      rev_diag_index_ptr_ = nullptr;
    }
  }
  if (!index_files_exist) {
    std::vector<std::vector<double>> shift_2d =
        proteoform_util::getNTermShift2D(
            proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    std::vector<std::vector<double>> n_term_acet_2d =
        proteoform_util::getNTermAcet2D(
            proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    term_index_ptr_ = mass_match_factory::getPrmTermMassMatchPtr(
        proteo_ptrs, shift_2d, mng_ptr->max_proteoform_mass_,
        mng_ptr->filter_scale_);
    diag_index_ptr_ = mass_match_factory::getPrmDiagMassMatchPtr(
        proteo_ptrs, mng_ptr->max_proteoform_mass_, mng_ptr->filter_scale_);
    std::vector<std::vector<double>> rev_shift_2d;
    std::vector<double> shift_1d(1, 0);
    for (size_t i = 0; i < proteo_ptrs.size(); i++) {
      rev_shift_2d.push_back(shift_1d);
    }
    rev_term_index_ptr_ = mass_match_factory::getSrmTermMassMatchPtr(
        proteo_ptrs, rev_shift_2d, n_term_acet_2d,
        mng_ptr->max_proteoform_mass_, mng_ptr->filter_scale_);
    rev_diag_index_ptr_ = mass_match_factory::getSrmDiagMassMatchPtr(
        proteo_ptrs, n_term_acet_2d, mng_ptr->max_proteoform_mass_,
        mng_ptr->filter_scale_);
  };
  term_nonincr_shifts_ =
      mass_match_util::compNonincreasingShifts(term_index_ptr_);
  diag_nonincr_shifts_ =
      mass_match_util::compNonincreasingShifts(diag_index_ptr_);
  std::vector<int> proteo_lens;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    proteo_lens.push_back(proteo_ptrs[i]->getResSeqPtr()->getLen());
  }
  row_ranges_.init(proteo_lens, term_index_ptr_, diag_index_ptr_,
                   rev_term_index_ptr_, rev_diag_index_ptr_,
                   prsm_para_ptr->getProtModPtrVec());
}

FilterSpecPtr OnePtmFilter::geneFilterSpec(
    const PrmMsPtrVec& prm_ms_ptr_vec, const PrmMsPtrVec& srm_ms_ptr_vec,
    const OnePtmFilterMngPtr& mng_ptr) {
  PeakTolerancePtr tole_ptr =
      mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  std::vector<std::pair<int, int>> pref_mass_errors =
      prm_ms_util::getIntMassErrorList(prm_ms_ptr_vec, tole_ptr,
                                       mng_ptr->filter_scale_, true, false);
  std::vector<std::pair<int, int>> suff_mass_errors =
      prm_ms_util::getIntMassErrorList(srm_ms_ptr_vec, tole_ptr,
                                       mng_ptr->filter_scale_, false, true);
  return std::make_shared<FilterSpec>(
      prm_ms_ptr_vec[0]->getMsHeaderPtr(),
      static_cast<int>(prm_ms_ptr_vec.size()), std::move(pref_mass_errors),
      std::move(suff_mass_errors), std::pair<int, int>(0, 0));
}

void OnePtmFilter::computeBestMatch(const PrmMsPtrVec& prm_ms_ptr_vec,
                                    const PrmMsPtrVec& srm_ms_ptr_vec) {
  computeBestMatch(geneFilterSpec(prm_ms_ptr_vec, srm_ms_ptr_vec, mng_ptr_));
}

void OnePtmFilter::computeBestMatch(const FilterSpecPtr& spec_ptr) {
  PeakTolerancePtr tole_ptr =
      mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  const std::vector<std::pair<int, int>>& pref_mass_errors =
      spec_ptr->getPrefMassErrors();
  const std::vector<std::pair<int, int>>& suff_mass_errors =
      spec_ptr->getSuffMassErrors();

  term_scores_.assign(term_index_ptr_->getRowNum(), 0);
  term_index_ptr_->compScores(pref_mass_errors, term_scores_);

  diag_scores_.assign(diag_index_ptr_->getRowNum(), 0);
  diag_index_ptr_->compScores(pref_mass_errors, diag_scores_);

  rev_term_scores_.assign(rev_term_index_ptr_->getRowNum(), 0);
  rev_term_index_ptr_->compScores(suff_mass_errors, rev_term_scores_);

  rev_diag_scores_.assign(rev_diag_index_ptr_->getRowNum(), 0);
  rev_diag_index_ptr_->compScores(suff_mass_errors, rev_diag_scores_);

  // the best row score of each proteoform in each index, computed once and
  // shared by the four proteoform types below
  mass_match_util::compProteoMaxScores(term_scores_, term_index_ptr_,
                                       term_max_);
  mass_match_util::compProteoMaxScores(diag_scores_, diag_index_ptr_,
                                       diag_max_);
  mass_match_util::compProteoMaxScores(rev_term_scores_, rev_term_index_ptr_,
                                       rev_term_max_);
  mass_match_util::compProteoMaxScores(rev_diag_scores_, rev_diag_index_ptr_,
                                       rev_diag_max_);

  double prec_minus_water_mass =
      spec_ptr->getHeaderPtr()->getFirstPrecMonoMassMinusWater();
  double prec_error_tole = tole_ptr->compPrecErrorTole(
      spec_ptr->getHeaderPtr()->getFirstPrecMonoMass());
  int group_spec_num = spec_ptr->getGroupSpecNum();

  // The candidates of a proteoform type (issue #348): with
  // range_candidates_, the best ranges, a protein once for each of its
  // ranges, each with the two shifts of its range; otherwise the best
  // proteins, each with the range of its best pair and the shifts of its
  // best pairs.
  auto findCands = [&](std::vector<short>& scores,
                       std::vector<short>& rev_scores,
                       const std::vector<int>& max_scores,
                       const std::vector<int>& rev_max_scores,
                       const MassMatchPtr& index_ptr,
                       const MassMatchPtr& rev_index_ptr, bool term_start,
                       bool term_end, int shift_num, int top_num,
                       const std::vector<char>* nonincr_shifts) {
    mass_match_util::RangeRows range_rows =
        row_ranges_.getRangeRows(term_start, term_end, false, false);
    if (mng_ptr_->range_candidates_) {
      return mass_match_util::findOneShiftTopRanges(
          scores, rev_scores, max_scores, rev_max_scores, index_ptr,
          rev_index_ptr, range_rows, prec_minus_water_mass, prec_error_tole,
          mng_ptr_->min_shift_, mng_ptr_->max_shift_, shift_num,
          mng_ptr_->threshold_, top_num, nonincr_shifts);
    }
    return mass_match_util::findOneShiftTopProteins(
        scores, rev_scores, max_scores, rev_max_scores, index_ptr,
        rev_index_ptr, prec_minus_water_mass, prec_error_tole,
        mng_ptr_->min_shift_, mng_ptr_->max_shift_, shift_num,
        mng_ptr_->threshold_, top_num, nonincr_shifts, &range_rows);
  };
  // the simple PrSM of a candidate, with the shift lists its type uses
  auto geneSimplePrsm = [&](const ProtCandidatePtr& cand_ptr, bool n_shifts,
                            bool c_shifts) {
    int id = cand_ptr->getProteinId();
    SimplePrsmPtr prsm_ptr = std::make_shared<SimplePrsm>(
        spec_ptr->getHeaderPtr(), group_spec_num, proteo_ptrs_[id],
        cand_ptr->getScore());
    prsm_ptr->setSeqRange(cand_ptr->getStartPos(), cand_ptr->getEndPos());
    if (n_shifts) {
      prsm_ptr->setNTermShifts(cand_ptr->getNTermShifts());
    }
    if (c_shifts) {
      prsm_ptr->setCTermShifts(cand_ptr->getCTermShifts());
    }
    return prsm_ptr;
  };

  ProtCandidatePtrVec comp_prots =
      findCands(term_scores_, rev_term_scores_, term_max_, rev_term_max_,
                term_index_ptr_, rev_term_index_ptr_, true, true,
                mng_ptr_->comp_cand_shift_num_, mng_ptr_->comp_num_,
                &term_nonincr_shifts_);
  comp_match_ptrs_.clear();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    comp_match_ptrs_.push_back(geneSimplePrsm(comp_prots[i], false, false));
  }

  ProtCandidatePtrVec pref_prots =
      findCands(term_scores_, rev_diag_scores_, term_max_, rev_diag_max_,
                term_index_ptr_, rev_diag_index_ptr_, true, false,
                mng_ptr_->cand_shift_num_, mng_ptr_->pref_suff_num_,
                &term_nonincr_shifts_);
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_prots.size(); i++) {
    pref_match_ptrs_.push_back(geneSimplePrsm(pref_prots[i], false, true));
  }

  ProtCandidatePtrVec suff_prots =
      findCands(diag_scores_, rev_term_scores_, diag_max_, rev_term_max_,
                diag_index_ptr_, rev_term_index_ptr_, false, true,
                mng_ptr_->cand_shift_num_, mng_ptr_->pref_suff_num_,
                &diag_nonincr_shifts_);
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_prots.size(); i++) {
    suff_match_ptrs_.push_back(geneSimplePrsm(suff_prots[i], true, false));
  }

  ProtCandidatePtrVec internal_prots =
      findCands(diag_scores_, rev_diag_scores_, diag_max_, rev_diag_max_,
                diag_index_ptr_, rev_diag_index_ptr_, false, false,
                mng_ptr_->cand_shift_num_, mng_ptr_->inte_num_,
                &diag_nonincr_shifts_);
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    internal_match_ptrs_.push_back(
        geneSimplePrsm(internal_prots[i], true, true));
  }
}

} /* namespace toppic */
