#include "ptmsearch/ptm_slow_match.hpp"

namespace prot {

PtmSlowMatch::PtmSlowMatch(ProteoformPtr proteo_ptr, 
                           SpectrumSetPtr spectrum_set_ptr,
                           CompShiftLowMemPtr comp_shift_ptr, 
                           PtmMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  deconv_ms_ptr_ = spectrum_set_ptr->getDeconvMsPtr();
  ms_six_ptr_=spectrum_set_ptr->getMsSixPtr();
  ms_three_ptr_ = spectrum_set_ptr->getMsThreePtr();
  proteo_ptr_ = proteo_ptr;
  initPsAlign(comp_shift_ptr);
}

// initialize ps_align
inline void PtmSlowMatch::initPsAlign(CompShiftLowMemPtr comp_shift_ptr){
  double scale = mng_ptr_->ptm_fast_filter_scale_;
  // n term strict c term nonstrict
  std::pair<std::vector<int>, std::vector<int>> sp_masses_toles 
      = getIntMassErrorList(ms_six_ptr_,scale,true,false);

  std::vector<double> best_shifts = comp_shift_ptr->findBestShift(
      sp_masses_toles.first,
      sp_masses_toles.second,
      proteo_ptr_->getBpSpecPtr()->getScaledPrmMasses(scale),
      mng_ptr_->n_top_diagonals_,
      mng_ptr_->min_diagonal_gap_,
      scale);
  DiagonalHeaderPtrVec n_term_shift_header_ptrs 
      = getNTermShiftHeaders (best_shifts, ms_six_ptr_, proteo_ptr_, mng_ptr_);

  /*
  for (size_t i = 0; i < n_term_shift_header_ptrs.size(); i++) {
    std::cout << " diagonal " << i << " shift " << n_term_shift_header_ptrs[i]->getProtNTermShift() << std::endl;
  }
  */

  BasicDiagonalPtrVec diagonal_ptrs = getDiagonals(n_term_shift_header_ptrs,
                                                   ms_six_ptr_,proteo_ptr_,mng_ptr_);
  for (size_t i = 0; i < diagonal_ptrs.size(); i++) {
    diagonal_ptrs[i]->getHeader()->setId(i);
  }
  std::vector<double> ms_masses = getMassList(ms_six_ptr_);
  std::vector<double> seq_masses = proteo_ptr_->getBpSpecPtr()->getPrmMasses();
  
  ps_align_ptr_ = PSAlignPtr(new PSAlign(ms_masses,seq_masses,diagonal_ptrs,mng_ptr_));
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

// get the header corresponding to the top left corner in the spectral grid
inline DiagonalHeaderPtr getTopLeftCornerHeader() {
  double shift = 0;
  // n_term strict; c_term nostrict; prot n_term match; prot c_term no_match
  // pep n_term no_match; pep c_term no_match
  return DiagonalHeaderPtr(
      new DiagonalHeader(shift, true, false, true, false, false, false));
}

inline DiagonalHeaderPtr getBottomRightCornerHeader(ProteoformPtr proteo_ptr,
                                                    PrmMsPtr ms_six_ptr) {
  double prec_mass = ms_six_ptr->getHeaderPtr()->getPrecMonoMass();
  double mole_mass = proteo_ptr->getResSeqPtr()->getSeqMass();
  double shift = prec_mass - mole_mass;
  // n term nostrict, c_term strict, prot n_term no match ; prot c_term match
  // pep n_term no match, pep c_term no match 
  return DiagonalHeaderPtr(
      new DiagonalHeader(shift, false, true, false, true, false, false));
}

inline int findSimilarShiftPos(const std::vector<double> &shifts, double s) {
  int best_pos = -1;
  double best_diff = std::numeric_limits<double>::infinity();
  for(size_t i = 0; i < shifts.size();i++){
    if(std::abs(shifts[i] - s) < best_diff){
      best_pos = i;
      best_diff = std::abs(shifts[i] - s);
    }
  }
  return best_pos;
}

inline bool isExist(const DiagonalHeaderPtrVec &header_ptrs, double shift) {
  for(size_t i = 0; i < header_ptrs.size();i++){
    if(header_ptrs[i]->getProtNTermShift() == shift) {
      return true;
    }
  }
  return false;
}

inline DiagonalHeaderPtrVec PtmSlowMatch::getNTermShiftHeaders(
    const std::vector<double> &best_shifts,
    PrmMsPtr ms_six_ptr,
    ProteoformPtr proteo_ptr,
    PtmMngPtr mng_ptr){

  DiagonalHeaderPtrVec header_ptrs = getNTermShiftListCommonHeaders(best_shifts);

  DiagonalHeaderPtrVec n_extend_header_ptrs;
  // get top-left corner header in spectral grid (shift is 0)
  DiagonalHeaderPtr top_left_corner_header_ptr = getTopLeftCornerHeader();
  n_extend_header_ptrs.push_back(top_left_corner_header_ptr);

  DiagonalHeaderPtrVec c_extend_header_ptrs;
  // get bottom-right corner header in the spectral grid. 
  DiagonalHeaderPtr bottom_right_corner_header_ptr 
      = getBottomRightCornerHeader(proteo_ptr, ms_six_ptr);
  c_extend_header_ptrs.push_back(bottom_right_corner_header_ptr);

  double prec_mass_minus_water = ms_six_ptr->getHeaderPtr()->getPrecMonoMass() - MassConstant::getWaterMass();
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
      if (!isExist(n_extend_header_ptrs, new_shift)) {
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
      if (!isExist(c_extend_header_ptrs, new_shift)) {
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

  return header_ptrs;
}


void PtmSlowMatch::compute(SemiAlignTypePtr type_ptr, PrsmPtrVec &prsm_ptrs) {
  ps_align_ptr_->compute(type_ptr);
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
    PrsmPtr prsm_ptr = geneResult(s);
    prsm_ptrs.push_back(prsm_ptr);
  }
}


PrsmPtr PtmSlowMatch::geneResult(int shift_num){
  DiagonalHeaderPtrVec header_ptrs= ps_align_ptr_->getResult(shift_num);
  //double score = ps_align_ptr_->getAlignScr(shift_num);
  //LOG_DEBUG("Shift " << shift_num << " score " << score);
  if(header_ptrs.size()==0) {
    return nullptr;
  }
  int first_pos = header_ptrs[0]->getTruncFirstResPos();
  int last_pos = header_ptrs[header_ptrs.size()-1]->getTruncLastResPos();
  ProteoformPtr sub_proteo_ptr  = getSubProteoform(proteo_ptr_, first_pos, last_pos);
  double refine_prec_mass = refinePrecursorAndHeaderShift(proteo_ptr_, ms_three_ptr_, 
                                                          header_ptrs, mng_ptr_);
  //double refine_prec_mass = ms_three_ptr_->getHeaderPtr()->getPrecMonoMass();

  double delta = refine_prec_mass - deconv_ms_ptr_->getHeaderPtr()->getPrecMonoMass();
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  ExtendMsPtr refine_ms_ptr_ = createMsThreePtr(deconv_ms_ptr_, delta, sp_para_ptr);

  DiagonalHeaderPtrVec refined_header_ptrs = refineHeadersBgnEnd(
      proteo_ptr_, refine_ms_ptr_, header_ptrs, mng_ptr_);

  if(refined_header_ptrs.size()==0){
    return nullptr;
  }

  ChangePtrVec changes = getUnexpectedChanges(refined_header_ptrs, 
                                              first_pos, last_pos);
  sub_proteo_ptr->addUnexpectedChangePtrVec(changes);

  return PrsmPtr(new Prsm(sub_proteo_ptr, deconv_ms_ptr_, refine_prec_mass,
          0, mng_ptr_->prsm_para_ptr_->getSpParaPtr()));
}

} /* namespace prot */
