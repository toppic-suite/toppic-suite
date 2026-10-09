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

#include "common/base/amino_acid.hpp"
#include "common/util/file_util.hpp"
#include "common/util/stage_profile.hpp"
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

// The rows of the diagonal index whose first residue can carry a truncation
// N-terminal modification (a variable PTM of the proteoform N-terminus) with
// the residues res_ptrs: the row starts at the third residue of the protein
// or later (the first two are those of the protein N-terminal forms, which
// the N-terminal form index has), and the amino acid there is one of theirs
// (Proteoform::setTruncNTermModPtr).
std::vector<char> compNTermRowMask(const MassMatchPtr& index_ptr,
                                   const std::vector<int>& row_starts,
                                   const ProteoformPtrVec& proteo_ptrs,
                                   const ResiduePtrVec& res_ptrs) {
  const std::vector<int>& row_proteo_ids = index_ptr->getRowProteoIds();
  std::vector<char> mask(row_proteo_ids.size(), 0);
  for (size_t r = 0; r < mask.size(); r++) {
    const ResSeqPtr& seq_ptr = proteo_ptrs[row_proteo_ids[r]]->getResSeqPtr();
    int start = row_starts[r];
    if (start < 2 || start >= seq_ptr->getLen()) {
      continue;
    }
    const AminoAcidPtr& acid_ptr =
        seq_ptr->getResiduePtr(start)->getAminoAcidPtr();
    for (size_t i = 0; i < res_ptrs.size(); i++) {
      if (res_ptrs[i]->getAminoAcidPtr() == acid_ptr) {
        mask[r] = 1;
        break;
      }
    }
  }
  return mask;
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
      stage_profile::Timer load_timer("var.load");
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
  initRows();
}

VarPtmFilter::VarPtmFilter(const ProteoformPtrVec& proteo_ptrs,
                           const VarPtmFilterMngPtr& mng_ptr,
                           MassMatchPtr term_index_ptr,
                           MassMatchPtr diag_index_ptr,
                           MassMatchPtr rev_term_index_ptr,
                           MassMatchPtr rev_diag_index_ptr)
    : mng_ptr_(mng_ptr),
      proteo_ptrs_(proteo_ptrs),
      diag_index_ptr_(std::move(diag_index_ptr)),
      rev_diag_index_ptr_(std::move(rev_diag_index_ptr)),
      term_index_ptr_(std::move(term_index_ptr)),
      rev_term_index_ptr_(std::move(rev_term_index_ptr)) {
  initRows();
}

void VarPtmFilter::initRows() {
  std::vector<int> proteo_lens;
  for (size_t i = 0; i < proteo_ptrs_.size(); i++) {
    proteo_lens.push_back(proteo_ptrs_[i]->getResSeqPtr()->getLen());
  }
  row_ranges_.init(proteo_lens, term_index_ptr_, diag_index_ptr_,
                   rev_term_index_ptr_, rev_diag_index_ptr_,
                   mng_ptr_->prsm_para_ptr_->getProtModPtrVec());
  for (size_t k = 0; k < mng_ptr_->n_term_shift_list_.size(); k++) {
    n_term_diag_masks_.push_back(
        compNTermRowMask(diag_index_ptr_, row_ranges_.diag_starts_,
                         proteo_ptrs_, mng_ptr_->n_term_res_ptr_vec_2d_[k]));
  }
}

FilterSpecPtr VarPtmFilter::geneFilterSpec(
    const PrmMsPtrVec& prm_ms_ptr_vec, const PrmMsPtrVec& srm_ms_ptr_vec,
    const VarPtmFilterMngPtr& mng_ptr) {
  const SpParaPtr& sp_para_ptr = mng_ptr->prsm_para_ptr_->getSpParaPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPrecTolerancePtr();
  SpPara::FilterMassListMode mode = sp_para_ptr->getFilterMassListMode();
  std::vector<std::pair<int, int>> pref_mass_errors =
      prm_ms_util::getIntMassErrorList(prm_ms_ptr_vec, tole_ptr,
                                       mng_ptr->filter_scale_, true, false,
                                       mode);
  std::vector<std::pair<int, int>> suff_mass_errors =
      prm_ms_util::getIntMassErrorList(srm_ms_ptr_vec, tole_ptr,
                                       mng_ptr->filter_scale_, false, true,
                                       mode);
  MsHeaderPtr header_ptr = prm_ms_ptr_vec[0]->getMsHeaderPtr();
  std::pair<int, int> prec_minus_water_mass_error =
      header_ptr->getFirstPrecMonoMassMinusWaterError(
          tole_ptr->compPrecErrorTole(header_ptr->getFirstPrecMonoMass()),
          mng_ptr->filter_scale_);
  return std::make_shared<FilterSpec>(
      prm_ms_ptr_vec[0]->getMsHeaderPtr(),
      static_cast<int>(prm_ms_ptr_vec.size()), std::move(pref_mass_errors),
      std::move(suff_mass_errors), prec_minus_water_mass_error);
}

void VarPtmFilter::scoreIndexes(
    const FilterSpecPtr& spec_ptr,
    const std::vector<std::pair<int, int>>& prec_minus_water_mass_errors,
    bool score_term, const std::function<void()>& base_scores_hook) {
  stage_profile::Timer score_timer("var.score");
  const std::vector<std::pair<int, int>>& pref_mass_errors =
      spec_ptr->getPrefMassErrors();
  const std::vector<std::pair<int, int>>& suff_mass_errors =
      spec_ptr->getSuffMassErrors();

  // the mass list scores (compMatchScores = compScores + compPrecScores)
  if (score_term) {
    term_scores_.assign(term_index_ptr_->getRowNum(), 0);
    term_index_ptr_->compScores(pref_mass_errors, term_scores_);
  }
  diag_scores_.assign(diag_index_ptr_->getRowNum(), 0);
  diag_index_ptr_->compScores(pref_mass_errors, diag_scores_);
  rev_term_scores_.assign(rev_term_index_ptr_->getRowNum(), 0);
  rev_term_index_ptr_->compScores(suff_mass_errors, rev_term_scores_);
  rev_diag_scores_.assign(rev_diag_index_ptr_->getRowNum(), 0);
  rev_diag_index_ptr_->compScores(suff_mass_errors, rev_diag_scores_);
  if (base_scores_hook) {
    base_scores_hook();
  }

  // the precursor bonus
  if (score_term) {
    term_prec_rows_.clear();
    term_index_ptr_->compPrecScores(prec_minus_water_mass_errors,
                                    term_scores_, term_prec_rows_);
  }
  diag_prec_rows_.clear();
  diag_index_ptr_->compPrecScores(prec_minus_water_mass_errors, diag_scores_,
                                  diag_prec_rows_);
  rev_term_prec_rows_.clear();
  rev_term_index_ptr_->compPrecScores(prec_minus_water_mass_errors,
                                      rev_term_scores_, rev_term_prec_rows_);
  rev_diag_prec_rows_.clear();
  rev_diag_index_ptr_->compPrecScores(prec_minus_water_mass_errors,
                                      rev_diag_scores_, rev_diag_prec_rows_);

  // the bonus rows of each proteoform, collected once and shared by the four
  // proteoform types
  if (score_term) {
    mass_match_util::compPrecMatchRows(term_index_ptr_, term_prec_rows_,
                                       term_rows_2d_);
  }
  mass_match_util::compPrecMatchRows(diag_index_ptr_, diag_prec_rows_,
                                     diag_rows_2d_);
  mass_match_util::compPrecMatchRows(rev_term_index_ptr_, rev_term_prec_rows_,
                                     rev_term_rows_2d_);
  mass_match_util::compPrecMatchRows(rev_diag_index_ptr_, rev_diag_prec_rows_,
                                     rev_diag_rows_2d_);
}

std::array<ProtCandidatePtrVec, 4> VarPtmFilter::findTypeCands(
    const FilterSpecPtr& spec_ptr, std::vector<double>& ptm_shifts,
    const std::vector<char>* diag_mask) {
  stage_profile::Timer find_timer("var.find");
  PeakTolerancePtr tole_ptr =
      mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  double prec_minus_water_mass =
      spec_ptr->getHeaderPtr()->getFirstPrecMonoMassMinusWater();
  double prec_error_tole = tole_ptr->compPrecErrorTole(
      spec_ptr->getHeaderPtr()->getFirstPrecMonoMass());

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
    const std::vector<char>* start_row_mask = term_start ? nullptr : diag_mask;
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
  // a truncation N-terminal modification is on a proteoform without a
  // protein N-terminal form, so its round has no complete and prefix
  // candidates
  if (diag_mask == nullptr) {
    cands[0] = findCands(term_scores_, rev_term_scores_, term_rows_2d_,
                         rev_term_rows_2d_, term_index_ptr_,
                         rev_term_index_ptr_, true, true,
                         mng_ptr_->comp_threshold_, mng_ptr_->comp_num_);
    cands[1] = findCands(term_scores_, rev_diag_scores_, term_rows_2d_,
                         rev_diag_rows_2d_, term_index_ptr_,
                         rev_diag_index_ptr_, true, false,
                         mng_ptr_->pref_suff_threshold_,
                         mng_ptr_->pref_suff_num_);
  }
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

void VarPtmFilter::computeBestMatch(
    const FilterSpecPtr& spec_ptr,
    const std::function<void()>& base_scores_hook) {
  scoreIndexes(spec_ptr,
               getShiftedMassErrors(spec_ptr->getPrecMinusWaterMassError(),
                                    mng_ptr_->int_shift_list_, false),
               true, base_scores_hook);
  std::array<ProtCandidatePtrVec, 4> cands =
      findTypeCands(spec_ptr, mng_ptr_->shift_list_, nullptr);

  int group_spec_num = spec_ptr->getGroupSpecNum();
  auto geneMatchPtrs = [&](const ProtCandidatePtrVec& prots,
                           SimplePrsmPtrVec* match_ptrs) {
    match_ptrs->clear();
    for (size_t i = 0; i < prots.size(); i++) {
      int id = prots[i]->getProteinId();
      match_ptrs->push_back(std::make_shared<SimplePrsm>(
          spec_ptr->getHeaderPtr(), group_spec_num, proteo_ptrs_[id],
          prots[i]));
    }
  };
  geneMatchPtrs(cands[0], &comp_match_ptrs_);
  geneMatchPtrs(cands[1], &pref_match_ptrs_);
  geneMatchPtrs(cands[2], &suff_match_ptrs_);
  geneMatchPtrs(cands[3], &internal_match_ptrs_);

  // One more round for each truncation N-terminal modification mass m (the
  // variable PTMs of the proteoform N-terminus): the spectrum with m taken
  // off its N-terminal masses and its precursor mass, so that the N-terminal
  // rows match the fragments of a proteoform whose N-terminus carries the
  // PTM; the other variable PTMs are those of up to var_ptm_num - 1, no
  // other PTM included. Only the diagonal index rows whose first residue
  // can carry the PTM start a pair (suffix and internal candidates). The
  // candidates hold the spectrum's own header (its precursor mass) and m,
  // the shift the search adds to them.
  const FilterSpecPtrVec& n_term_specs = spec_ptr->getNTermModSpecs();
  n_term_match_ptrs_.resize(n_term_specs.size());
  for (size_t k = 0; k < n_term_specs.size(); k++) {
    const FilterSpecPtr& n_term_spec_ptr = n_term_specs[k];
    scoreIndexes(n_term_spec_ptr,
                 getShiftedMassErrors(
                     n_term_spec_ptr->getPrecMinusWaterMassError(),
                     mng_ptr_->n_term_other_int_shift_list_, true),
                 false);
    std::array<ProtCandidatePtrVec, 4> round_cands =
        findTypeCands(n_term_spec_ptr, mng_ptr_->n_term_other_shift_list_,
                      &n_term_diag_masks_[k]);
    std::vector<double> trunc_shifts(1, mng_ptr_->n_term_shift_list_[k]);
    for (size_t t = 0; t < round_cands.size(); t++) {
      geneMatchPtrs(round_cands[t], &n_term_match_ptrs_[k][t]);
      for (size_t i = 0; i < n_term_match_ptrs_[k][t].size(); i++) {
        n_term_match_ptrs_[k][t][i]->setTruncNTermModShifts(trunc_shifts);
      }
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
