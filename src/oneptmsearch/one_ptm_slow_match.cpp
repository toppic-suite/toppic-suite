#include "ptmsearch/diagonal_header_util.hpp"
#include "oneptmsearch/one_ptm_slow_match.hpp"

namespace prot {

OnePtmSlowMatch::OnePtmSlowMatch(ProteoformPtr proteo_ptr, 
                           SpectrumSetPtr spectrum_set_ptr,
                           CompShiftLowMemPtr comp_shift_ptr, 
                           SemiAlignTypePtr align_type_ptr,
                           OnePtmSearchMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  align_type_ptr_ = align_type_ptr;
  prec_mono_mass_ = spectrum_set_ptr->getPrecMonoMass();
  deconv_ms_ptr_vec_ = spectrum_set_ptr->getDeconvMsPtrVec();
  ms_six_ptr_vec_ = spectrum_set_ptr->getMsSixPtrVec();
  ms_three_ptr_vec_ = spectrum_set_ptr->getMsThreePtrVec();

  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  prm_peaks_ = getPrmPeakPtrs(ms_six_ptr_vec_, tole_ptr);
  group_spec_num_ = ms_six_ptr_vec_.size();

  proteo_ptr_ = proteo_ptr;
  initOnePtmAlign(comp_shift_ptr);
}

// initialize ps_align
inline void OnePtmSlowMatch::initOnePtmAlign(CompShiftLowMemPtr comp_shift_ptr){
  BasicDiagonalPtrVec diagonal_ptrs = compAlignDiagonals(comp_shift_ptr);
  std::vector<double> seq_masses = proteo_ptr_->getBpSpecPtr()->getPrmMasses();
  std::vector<double> ms_masses;
  for (size_t i = 0; i < prm_peaks_.size(); i++) {
    ms_masses.push_back(prm_peaks_[i]->getPosition());
  }
  ps_align_ptr_ = PSAlignPtr(new PSAlign(ms_masses,seq_masses,diagonal_ptrs,mng_ptr_->align_mng_ptr_));
}

inline std::vector<double> OnePtmSlowMatch::compBestShifts(CompShiftLowMemPtr comp_shift_ptr) {
  double scale = mng_ptr_->ptm_fast_filter_scale_;
  // n term strict c term nonstrict
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> sp_masses_toles 
      = getIntMassErrorList(ms_six_ptr_vec_,tole_ptr,scale,true,false);
  return comp_shift_ptr->findBestShift(
      sp_masses_toles,
      proteo_ptr_->getBpSpecPtr()->getScaledPrmMasses(scale),
      mng_ptr_->n_top_diagonals_,
      mng_ptr_->min_diagonal_gap_,
      scale);
}

inline void OnePtmSlowMatch::extendNHeaders(std::vector<double> &best_shifts,
                                            DiagonalHeaderPtrVec &n_extend_header_ptrs) {
  std::vector<double> prm_masses = proteo_ptr_->getBpSpecPtr()->getPrmMasses();
  // shifts for n_term matches
  std::vector<double> n_term_match_shifts;
  for(size_t i=1;i< prm_masses.size();i++){
    n_term_match_shifts.push_back(-prm_masses[i]);
  }
  // add trunc headers that have similar shift to best shift headers
  for (size_t i = 0; i < best_shifts.size(); i++) {
    double s = best_shifts[i];
    // find a similar shift in n_term_match_shifts
    int best_n_pos = findSimilarShiftPos(n_term_match_shifts, s);
    if (best_n_pos >= 0) {
      double new_shift = n_term_match_shifts[best_n_pos];
      if (!isExistHeader(n_extend_header_ptrs, new_shift)) {
        // n_term strict; c_term nostrict; prot n_term no_match; prot c_term no_match
        // pep n_term match; pep c_term no_match
        DiagonalHeaderPtr header_ptr(
            new DiagonalHeader(new_shift, true, false, false, false, true, false));
        n_extend_header_ptrs.push_back(header_ptr);
      }
    }
  }
}

inline void OnePtmSlowMatch::extendCHeaders(std::vector<double> &best_shifts,
                                                  DiagonalHeaderPtrVec &c_extend_header_ptrs) {
  std::vector<double> prm_masses = proteo_ptr_->getBpSpecPtr()->getPrmMasses();
  std::vector<double> c_term_match_shifts;
  double prec_mass_minus_water = prec_mono_mass_ - MassConstant::getWaterMass();
  for(size_t i=1;i< prm_masses.size();i++){
    c_term_match_shifts.push_back(prec_mass_minus_water - prm_masses[i]);
  }
  for (size_t i = 0; i < best_shifts.size(); i++) {
    double s = best_shifts[i];
    int best_c_pos = findSimilarShiftPos(c_term_match_shifts, s);
    //LOG_DEBUG("Shift " << s <<" C term position " << best_c_pos);
    if (best_c_pos >= 0) {
      double new_shift = c_term_match_shifts[best_c_pos];
      //LOG_DEBUG("Shift " << s <<" C term shift " << new_shift);
      if (!isExistHeader(c_extend_header_ptrs, new_shift)) {
        // n term nostrict, c_term strict, prot n_term no match ; prot c_term no match
        // pep n_term no match, pep c_term match 
        DiagonalHeaderPtr header_ptr(
            new DiagonalHeader(new_shift, false, true, false, false, false, true));
        c_extend_header_ptrs.push_back(header_ptr);
      }
    }
  }
}

inline BasicDiagonalPtrVec OnePtmSlowMatch::getDiagonals(DiagonalHeaderPtrVec &header_ptrs) {
  BasicDiagonalPtrVec results;
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    double n_shift = header_ptrs[i]->getProtNTermShift();
    double c_shift = prec_mono_mass_ - proteo_ptr_->getResSeqPtr()->getSeqMass() - n_shift;
    header_ptrs[i]->initHeader(c_shift, proteo_ptr_, mng_ptr_->align_prefix_suffix_shift_thresh_);
    BasicDiagonalPtrVec diagonal_ptrs = geneDiagonals(header_ptrs,
                                                      prm_peaks_, group_spec_num_,
                                                      proteo_ptr_);
    for (size_t j = 0; j < diagonal_ptrs.size(); j++) {
      if ((int)j < mng_ptr_->top_diag_num_ || (int)diagonal_ptrs[i]->size() >= mng_ptr_->diag_min_size_) {
        results.push_back(diagonal_ptrs[j]);
      }
    }
  }
  return results;
}

inline BasicDiagonalPtrVec OnePtmSlowMatch::compAlignDiagonals(CompShiftLowMemPtr comp_shift_ptr) {
  DiagonalHeaderPtrVec n_extend_header_ptrs;
  // get top-left corner header in spectral grid (shift is 0)
  DiagonalHeaderPtr top_left_corner_header_ptr = getTopLeftCornerHeader();
  n_extend_header_ptrs.push_back(top_left_corner_header_ptr);

  DiagonalHeaderPtrVec c_extend_header_ptrs;
  // get bottom-right corner header in the spectral grid. 
  double seq_mass = proteo_ptr_->getResSeqPtr()->getSeqMass();
  DiagonalHeaderPtr bottom_right_corner_header_ptr 
      = getBottomRightCornerHeader(seq_mass, prec_mono_mass_);
  c_extend_header_ptrs.push_back(bottom_right_corner_header_ptr);

  // get best shifts
  std::vector<double> best_shifts;
  if (align_type_ptr_ != SemiAlignTypeFactory::getCompletePtr()) {
    best_shifts = compBestShifts(comp_shift_ptr);
  }

  // add n terminal extended shifts 
  if (align_type_ptr_ == SemiAlignTypeFactory::getSuffixPtr() 
      || align_type_ptr_ == SemiAlignTypeFactory::getInternalPtr()) {
    extendNHeaders(best_shifts, n_extend_header_ptrs);
  }

  if (align_type_ptr_ == SemiAlignTypeFactory::getPrefixPtr() 
      || align_type_ptr_ == SemiAlignTypeFactory::getInternalPtr()) {
    extendCHeaders(best_shifts, c_extend_header_ptrs);
  }

  BasicDiagonalPtrVec diagonal_ptrs = getDiagonals(n_extend_header_ptrs);
  BasicDiagonalPtrVec c_diagonal_ptrs = getDiagonals(c_extend_header_ptrs);

  diagonal_ptrs.insert(diagonal_ptrs.end(), c_diagonal_ptrs.begin(), c_diagonal_ptrs.end());
  return diagonal_ptrs;
}

PrsmPtr OnePtmSlowMatch::getResultPrsm() {
  return ps_align_ptr_->geneResult(1, proteo_ptr_, deconv_ms_ptr_vec_, ms_three_ptr_vec_, mng_ptr_->prsm_para_ptr_);
}

} /* namespace prot */
