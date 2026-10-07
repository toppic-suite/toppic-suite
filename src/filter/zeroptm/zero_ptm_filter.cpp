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

#include "filter/zeroptm/zero_ptm_filter.hpp"

#include <atomic>
#include <filesystem>
#include <mutex>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/massmatch/index_load_gate.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "ms/factory/extend_ms_util.hpp"
#include "seq/prot_candidate.hpp"
#include "seq/proteoform_util.hpp"

namespace toppic {

// The first block of the stage is deserialized under the mutex, which
// initializes Boost.Serialization's lazily-constructed type-registration
// singletons on one thread; after that, concurrent binary_iarchive instances
// on distinct streams are safe, and the per-block loads run in parallel.
// Fully serializing every load, as this mutex used to, made the loads the
// stage's floor: on the 480 benchmark the 12 worker threads spent 823
// thread-seconds waiting on this lock around 98 s of one-at-a-time reading,
// so the 100 s stage was the serialized load channel itself (toppic_test
// experiment 66). Concurrency changes only when each block's vectors fill,
// never their contents.
namespace {
std::mutex zero_ptm_filter_mutex;
std::atomic<bool> zero_ptm_filter_warmed{false};
}  // namespace
ZeroPtmFilter::ZeroPtmFilter(const ProteoformPtrVec& proteo_ptrs,
                             const ZeroPtmFilterMngPtr& mng_ptr,
                             std::string block_str) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

  std::string index_dir = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx";

  std::string parameters = mng_ptr->getIndexFilePara();
  std::string suffix = parameters + block_str;  // last part of file name

  // check if all index files for this ptm is present. if not, generate index
  // files again.
  bool index_files_exist = true;
  for (size_t t = 0; t < mng_ptr->zero_ptm_file_vec_.size(); t++) {
    std::string file_name = mng_ptr->zero_ptm_file_vec_[t] + suffix;
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
      std::unique_lock<std::mutex> lock(zero_ptm_filter_mutex,
                                        std::defer_lock);
      if (!zero_ptm_filter_warmed.load(std::memory_order_acquire)) {
        lock.lock();
        // threads that queued behind the first load proceed in parallel
        if (zero_ptm_filter_warmed.load(std::memory_order_relaxed)) {
          lock.unlock();
        }
      }
      IndexLoadGate::Slot load_slot(index_load_gate);
      term_index_ptr_->deserializeMassMatch(
          mng_ptr->zero_ptm_file_vec_[0] + suffix, index_dir);
      diag_index_ptr_->deserializeMassMatch(
          mng_ptr->zero_ptm_file_vec_[1] + suffix, index_dir);
      rev_term_index_ptr_->deserializeMassMatch(
          mng_ptr->zero_ptm_file_vec_[2] + suffix, index_dir);
      rev_diag_index_ptr_->deserializeMassMatch(
          mng_ptr->zero_ptm_file_vec_[3] + suffix, index_dir);
      if (lock.owns_lock()) {
        zero_ptm_filter_warmed.store(true, std::memory_order_release);
      }
    }
  } else {
    LOG_DEBUG("get shifts");
    std::vector<std::vector<double> > shift_2d =
        proteoform_util::getNTermShift2D(
            proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    std::vector<std::vector<double> > n_term_acet_2d =
        proteoform_util::getNTermAcet2D(
            proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    LOG_DEBUG("get shifts complete");
    // N-terminal indexes
    term_index_ptr_ = mass_match_factory::getPrmTermMassMatchPtr(
        proteo_ptrs, shift_2d, mng_ptr->max_proteoform_mass_,
        mng_ptr->filter_scale_);
    // Prm indexes
    diag_index_ptr_ = mass_match_factory::getPrmDiagMassMatchPtr(
        proteo_ptrs, mng_ptr->max_proteoform_mass_, mng_ptr->filter_scale_);
    LOG_DEBUG("diag index");
    std::vector<std::vector<double> > rev_shift_2d;
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
                   prsm_para_ptr->getProtModPtrVec());
}

FilterSpecPtr ZeroPtmFilter::geneFilterSpec(
    const ExtendMsPtrVec& ms_ptr_vec, const ZeroPtmFilterMngPtr& mng_ptr) {
  PeakTolerancePtr tole_ptr =
      mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  bool pref = true;
  std::vector<std::pair<int, int> > pref_mass_errors =
      extend_ms_util::getExtendIntMassErrorList(ms_ptr_vec, pref,
                                                mng_ptr->filter_scale_);
  pref = false;
  std::vector<std::pair<int, int> > suff_mass_errors =
      extend_ms_util::getExtendIntMassErrorList(ms_ptr_vec, pref,
                                                mng_ptr->filter_scale_);
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  std::pair<int, int> prec_minus_water_mass_error =
      header_ptr->getFirstPrecMonoMassMinusWaterError(
          tole_ptr->compPrecErrorTole(header_ptr->getFirstPrecMonoMass()),
          mng_ptr->filter_scale_);
  return std::make_shared<FilterSpec>(
      ms_ptr_vec[0]->getMsHeaderPtr(), static_cast<int>(ms_ptr_vec.size()),
      std::move(pref_mass_errors), std::move(suff_mass_errors),
      prec_minus_water_mass_error);
}

void ZeroPtmFilter::computeBestMatch(const ExtendMsPtrVec& ms_ptr_vec) {
  computeBestMatch(geneFilterSpec(ms_ptr_vec, mng_ptr_));
}

void ZeroPtmFilter::computeBestMatch(const FilterSpecPtr& spec_ptr) {
  const std::vector<std::pair<int, int> >& pref_mass_errors =
      spec_ptr->getPrefMassErrors();
  const std::vector<std::pair<int, int> >& suff_mass_errors =
      spec_ptr->getSuffMassErrors();
  // one precursor mass; the rows that receive its bonus are noted
  const std::vector<std::pair<int, int> > prec_minus_water_mass_errors(
      1, spec_ptr->getPrecMinusWaterMassError());

  int term_row_num = term_index_ptr_->getRowNum();
  std::vector<short> term_scores(term_row_num, 0);
  term_prec_rows_.clear();
  term_index_ptr_->compMatchScores(pref_mass_errors,
                                   prec_minus_water_mass_errors, term_scores,
                                   term_prec_rows_);

  int diag_row_num = diag_index_ptr_->getRowNum();
  std::vector<short> diag_scores(diag_row_num, 0);
  diag_prec_rows_.clear();
  diag_index_ptr_->compMatchScores(pref_mass_errors,
                                   prec_minus_water_mass_errors, diag_scores,
                                   diag_prec_rows_);

  int rev_term_row_num = rev_term_index_ptr_->getRowNum();
  std::vector<short> rev_term_scores(rev_term_row_num, 0);
  rev_term_prec_rows_.clear();
  rev_term_index_ptr_->compMatchScores(suff_mass_errors,
                                       prec_minus_water_mass_errors,
                                       rev_term_scores, rev_term_prec_rows_);

  int rev_diag_row_num = rev_diag_index_ptr_->getRowNum();
  std::vector<short> rev_diag_scores(rev_diag_row_num, 0);
  rev_diag_prec_rows_.clear();
  rev_diag_index_ptr_->compMatchScores(suff_mass_errors,
                                       prec_minus_water_mass_errors,
                                       rev_diag_scores, rev_diag_prec_rows_);

  int threshold = MassMatch::getPrecursorMatchScore() * 2 + 4;
  int group_spec_num = spec_ptr->getGroupSpecNum();
  // The candidates are sequence ranges (issue #348): the pairs of rows that
  // both received the precursor bonus and whose range mass matches the
  // precursor mass, with the ranges the zero-shift search of the type takes
  // (PREFIX not to the last residue, SUFFIX not from the first, INTERNAL
  // neither). The bonus rows of each proteoform are collected once and
  // shared by the four proteoform types below.
  mass_match_util::compPrecMatchRows(term_index_ptr_, term_prec_rows_,
                                     term_rows_2d_);
  mass_match_util::compPrecMatchRows(diag_index_ptr_, diag_prec_rows_,
                                     diag_rows_2d_);
  mass_match_util::compPrecMatchRows(rev_term_index_ptr_, rev_term_prec_rows_,
                                     rev_term_rows_2d_);
  mass_match_util::compPrecMatchRows(rev_diag_index_ptr_, rev_diag_prec_rows_,
                                     rev_diag_rows_2d_);
  PeakTolerancePtr tole_ptr =
      mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  double prec_minus_water_mass =
      spec_ptr->getHeaderPtr()->getFirstPrecMonoMassMinusWater();
  double prec_error_tole = tole_ptr->compPrecErrorTole(
      spec_ptr->getHeaderPtr()->getFirstPrecMonoMass());
  bool one_per_protein = !mng_ptr_->range_candidates_;
  LOG_DEBUG("allow complete prot"
            << mng_ptr_->prsm_para_ptr_->allowCompleteProt());
  if (mng_ptr_->prsm_para_ptr_->allowCompleteProt()) {
    ProtCandidatePtrVec comp_prots = mass_match_util::findZeroShiftTopRanges(
        term_scores, rev_term_scores, term_rows_2d_, rev_term_rows_2d_,
        term_index_ptr_, rev_term_index_ptr_,
        row_ranges_.getRangeRows(true, true, false, false),
        prec_minus_water_mass, prec_error_tole, threshold, mng_ptr_->comp_num_,
        one_per_protein);
    comp_match_ptrs_.clear();
    for (size_t i = 0; i < comp_prots.size(); i++) {
      int id = comp_prots[i]->getProteinId();
      comp_match_ptrs_.push_back(std::make_shared<SimplePrsm>(
          spec_ptr->getHeaderPtr(), group_spec_num, proteo_ptrs_[id],
          comp_prots[i]));
    }
  }
  if (mng_ptr_->prsm_para_ptr_->allowPrefixProt()) {
    ProtCandidatePtrVec pref_prots = mass_match_util::findZeroShiftTopRanges(
        term_scores, rev_diag_scores, term_rows_2d_, rev_diag_rows_2d_,
        term_index_ptr_, rev_diag_index_ptr_,
        row_ranges_.getRangeRows(true, false, false, true),
        prec_minus_water_mass, prec_error_tole, threshold,
        mng_ptr_->pref_suff_num_, one_per_protein);
    pref_match_ptrs_.clear();
    for (size_t i = 0; i < pref_prots.size(); i++) {
      int id = pref_prots[i]->getProteinId();
      pref_match_ptrs_.push_back(std::make_shared<SimplePrsm>(
          spec_ptr->getHeaderPtr(), group_spec_num, proteo_ptrs_[id],
          pref_prots[i]));
    }
  }

  if (mng_ptr_->prsm_para_ptr_->allowSuffixProt()) {
    ProtCandidatePtrVec suff_prots = mass_match_util::findZeroShiftTopRanges(
        diag_scores, rev_term_scores, diag_rows_2d_, rev_term_rows_2d_,
        diag_index_ptr_, rev_term_index_ptr_,
        row_ranges_.getRangeRows(false, true, true, false),
        prec_minus_water_mass, prec_error_tole, threshold,
        mng_ptr_->pref_suff_num_, one_per_protein);
    suff_match_ptrs_.clear();
    for (size_t i = 0; i < suff_prots.size(); i++) {
      int id = suff_prots[i]->getProteinId();
      suff_match_ptrs_.push_back(std::make_shared<SimplePrsm>(
          spec_ptr->getHeaderPtr(), group_spec_num, proteo_ptrs_[id],
          suff_prots[i]));
    }
  }

  if (mng_ptr_->prsm_para_ptr_->allowInternalProt()) {
    ProtCandidatePtrVec internal_prots =
        mass_match_util::findZeroShiftTopRanges(
            diag_scores, rev_diag_scores, diag_rows_2d_, rev_diag_rows_2d_,
            diag_index_ptr_, rev_diag_index_ptr_,
            row_ranges_.getRangeRows(false, false, true, true),
            prec_minus_water_mass, prec_error_tole, threshold,
            mng_ptr_->inte_num_, one_per_protein);
    internal_match_ptrs_.clear();
    for (size_t i = 0; i < internal_prots.size(); i++) {
      int id = internal_prots[i]->getProteinId();
      internal_match_ptrs_.push_back(std::make_shared<SimplePrsm>(
          spec_ptr->getHeaderPtr(), group_spec_num, proteo_ptrs_[id],
          internal_prots[i]));
    }
  }
}

} /* namespace toppic */
