#include "ptmsearch/diagonal_header_util.hpp"
#include "ptmsearch/ptm_slow_match.hpp"

namespace prot {

PtmSlowMatch::PtmSlowMatch(ProteoformPtr proteo_ptr, 
                           SpectrumSetPtr spectrum_set_ptr,
                           CompShiftLowMemPtr comp_shift_ptr, 
                           PtmMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  prec_mono_mass_ = spectrum_set_ptr->getPrecMonoMass();
  deconv_ms_ptr_vec_ = spectrum_set_ptr->getDeconvMsPtrVec();
  ms_six_ptr_vec_ = spectrum_set_ptr->getMsSixPtrVec();
  ms_three_ptr_vec_ = spectrum_set_ptr->getMsThreePtrVec();
  proteo_ptr_ = proteo_ptr;
  initPsAlign(comp_shift_ptr);
}

// initialize ps_align
inline void PtmSlowMatch::initPsAlign(CompShiftLowMemPtr comp_shift_ptr){
  double scale = mng_ptr_->ptm_fast_filter_scale_;
  // n term strict c term nonstrict
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> sp_masses_toles 
      = PrmMs::getIntMassErrorList(ms_six_ptr_vec_,tole_ptr,scale,true,false);
  std::vector<double> best_shifts = comp_shift_ptr->findBestShift(
      sp_masses_toles,
      proteo_ptr_->getBpSpecPtr()->getScaledPrmMasses(scale),
      mng_ptr_->n_top_diagonals_,
      mng_ptr_->min_diagonal_gap_,
      scale);
  DiagonalHeaderPtrVec n_term_shift_header_ptrs 
      = getNTermShiftHeaders (best_shifts, prec_mono_mass_, proteo_ptr_, mng_ptr_);

  PrmPeakPtrVec prm_peaks = PrmMs::getPrmPeakPtrs(ms_six_ptr_vec_, tole_ptr);
  int group_spec_num = ms_six_ptr_vec_.size();
  BasicDiagonalPtrVec diagonal_ptrs = geneDiagonals(n_term_shift_header_ptrs,
                                                   prm_peaks, group_spec_num,
                                                   proteo_ptr_);
  std::vector<double> seq_masses = proteo_ptr_->getBpSpecPtr()->getPrmMasses();
  std::vector<double> ms_masses;
  for (size_t i = 0; i < prm_peaks.size(); i++) {
    ms_masses.push_back(prm_peaks[i]->getPosition());
  }

  ps_align_ptr_ = PSAlignPtr(new PSAlign(ms_masses, seq_masses,
                                         diagonal_ptrs, mng_ptr_->align_mng_ptr_));
}

// get headers without n trunc and c trunc information 
inline DiagonalHeaderPtrVec getNTermShiftListCommonHeaders(std::vector<double> best_shifts) {
  DiagonalHeaderPtrVec header_ptrs;
  for (size_t i = 0; i < best_shifts.size(); i++) {
    // n term shift; c term nostrict; no prot nterm match; no prot cterm match,
    // no pep nterm match; no pep cterm match
    DiagonalHeaderPtr header_ptr(new DiagonalHeader(best_shifts[i], 
                                                    true, false, false, false, false, false));
    header_ptrs.push_back(header_ptr);
  }
  return header_ptrs;
}

inline DiagonalHeaderPtrVec PtmSlowMatch::getNTermShiftHeaders(
    const std::vector<double> &best_shifts,
    double prec_mass,
    ProteoformPtr proteo_ptr,
    PtmMngPtr mng_ptr){

  DiagonalHeaderPtrVec header_ptrs = getNTermShiftListCommonHeaders(best_shifts);

  DiagonalHeaderPtrVec n_extend_header_ptrs;
  // get top-left corner header in spectral grid (shift is 0)
  DiagonalHeaderPtr top_left_corner_header_ptr = getTopLeftCornerHeader();
  n_extend_header_ptrs.push_back(top_left_corner_header_ptr);

  DiagonalHeaderPtrVec c_extend_header_ptrs;
  // get bottom-right corner header in the spectral grid. 
  double seq_mass = proteo_ptr_->getResSeqPtr()->getSeqMass();
  DiagonalHeaderPtr bottom_right_corner_header_ptr 
      = getBottomRightCornerHeader(seq_mass, prec_mass);
  c_extend_header_ptrs.push_back(bottom_right_corner_header_ptr);

  double prec_mass_minus_water = prec_mass - MassConstant::getWaterMass();
  std::vector<double> prm_masses = proteo_ptr_->getBpSpecPtr()->getPrmMasses();
  // shifts for n_term matches
  std::vector<double> n_term_match_shifts;
  // shifts for c_term matches
  std::vector<double> c_term_match_shifts;
  for(size_t i=1;i< prm_masses.size();i++){
    n_term_match_shifts.push_back(-prm_masses[i]);
    c_term_match_shifts.push_back(prec_mass_minus_water - prm_masses[i]);
  }

  // add trunc headers that have similar shift to best shift headers
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    double s = header_ptrs[i]->getProtNTermShift();
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
  header_ptrs.insert(header_ptrs.end(), n_extend_header_ptrs.begin(), n_extend_header_ptrs.end());
  header_ptrs.insert(header_ptrs.end(), c_extend_header_ptrs.begin(), c_extend_header_ptrs.end());

  for (size_t i = 0; i < header_ptrs.size(); i++) {
    double n_shift = header_ptrs[i]->getProtNTermShift();
    double c_shift = prec_mono_mass_ - proteo_ptr_->getResSeqPtr()->getSeqMass() - n_shift;
    header_ptrs[i]->initHeader(c_shift, proteo_ptr_, mng_ptr_->align_prefix_suffix_shift_thresh_);
  }

  return header_ptrs;
}

void PtmSlowMatch::compute(AlignTypePtr type_ptr, PrsmPtrVec &prsm_ptrs) {
    ps_align_ptr_->compute(type_ptr);
    for (int s = 1; s <= mng_ptr_->getShiftNum(); s++) {
        PrsmPtr prsm_ptr = ps_align_ptr_->geneResult(s, proteo_ptr_, deconv_ms_ptr_vec_, 
                ms_three_ptr_vec_, mng_ptr_->prsm_para_ptr_);
        prsm_ptrs.push_back(prsm_ptr);
    }
}


//PrsmPtr PtmSlowMatch::geneResult(int shift_num){
  //DiagonalHeaderPtrVec header_ptrs= ps_align_ptr_->getResult(shift_num);
  ////double score = ps_align_ptr_->getAlignScr(shift_num);
  ////LOG_DEBUG("Shift " << shift_num << " score " << score);
  //if(header_ptrs.size()==0) {
    //return nullptr;
  //}
  //int first_pos = header_ptrs[0]->getTruncFirstResPos();
  //int last_pos = header_ptrs[header_ptrs.size()-1]->getTruncLastResPos();
  //ProteoformPtr sub_proteo_ptr  = getSubProteoform(proteo_ptr_, first_pos, last_pos);

  //double refine_prec_mass = refinePrecursorAndHeaderShift(proteo_ptr_, ms_three_ptr_vec_, 
                                                          //header_ptrs, mng_ptr_);

  //SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  //ExtendMsPtrVec refine_ms_ptr_vec_ = createMsThreePtrVec(deconv_ms_ptr_vec_,  sp_para_ptr, refine_prec_mass);

  //DiagonalHeaderPtrVec refined_header_ptrs = refineHeadersBgnEnd(
      //proteo_ptr_, refine_ms_ptr_vec_, header_ptrs, mng_ptr_);

  //if(refined_header_ptrs.size()==0){
    //return nullptr;
  //}

  //ChangePtrVec changes = getUnexpectedChanges(refined_header_ptrs, 
                                              //first_pos, last_pos);
  //sub_proteo_ptr->addUnexpectedChangePtrVec(changes);

  //return PrsmPtr(new Prsm(sub_proteo_ptr, deconv_ms_ptr_vec_, refine_prec_mass,
          //mng_ptr_->prsm_para_ptr_->getSpParaPtr()));
/*}*/

} /* namespace prot */
