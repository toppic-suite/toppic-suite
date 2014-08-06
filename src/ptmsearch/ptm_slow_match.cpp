/*
 * ptm_slow_match.cpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#include "ptmsearch/ptm_slow_match.hpp"

namespace prot {

PtmSlowMatch::PtmSlowMatch(ProteoformPtr proteoform, SpectrumSetPtr spectrum_set,
                           CompShiftLowMemPtr comp_shift, PtmMngPtr mng){
  mng_=mng;
  deconv_ms_ = spectrum_set->getDeconvMsPtr();
  ms_six_=spectrum_set->getMsSixPtr();
  ms_three_ = spectrum_set->getMsThreePtr();
  proteoform_ = proteoform;
  initPsAlign(comp_shift);
}

// initialize ps_align
void PtmSlowMatch::initPsAlign(CompShiftLowMemPtr comp_shift){
  double scale = mng_->ptm_fast_filter_scale_;
  // n term strict c term nonstrict
  std::vector<std::vector<int>> sp_masses_toles = getIntMassErrorList(ms_six_,scale,true,false);

  std::vector<double> best_shifts = comp_shift->findBestShift(
      sp_masses_toles[0],
      sp_masses_toles[1],
      proteoform_->getBpSpecPtr()->getScaledPrmMasses(scale),
      mng_->n_top_diagonals_,
      mng_->min_diagonal_gap_,
      scale);
  DiagonalHeaderPtrVec n_term_shift_headers 
      = getNTermShiftHeaders (best_shifts, ms_six_, proteoform_, mng_);
  BasicDiagPairDiagPtrVec diagonals = getDiagonals(n_term_shift_headers,
                                                   ms_six_,proteoform_,mng_);
  /*
  for (unsigned int i = 0; i < diagonals.size(); i++) {
    std::cout << "diagonal " << i << " shift " << diagonals[i]->getHeader()->getProtNTermShift()
        << " match point num " << diagonals[i]->size() << std::endl;
  }
  */
  std::vector<double> ms_masses = getMassList(ms_six_);
  std::vector<double> seq_masses = proteoform_->getBpSpecPtr()->getPrmMasses();
  ps_align_ = PSAlignPtr(new PSAlign(ms_masses,seq_masses,diagonals,mng_));
}

// get headers without n trunc and c trunc information 
DiagonalHeaderPtrVec getNTermShiftListCommonHeaders(std::vector<double> best_shifts) {
  DiagonalHeaderPtrVec headers;
  for (unsigned int i = 0; i < best_shifts.size(); i++) {
    // n term shift; c term nostrict; no prot nterm match; no prot cterm match,
    // no pep nterm match; no pep cterm match
    DiagonalHeaderPtr header_ptr(new DiagonalHeader(best_shifts[i], 
                                                    true, false, false, false, false, false));
    headers.push_back(header_ptr);
  }
  return headers;
}

// get the header corresponding to the top left corner in the spectral grid
DiagonalHeaderPtr getTopLeftCornerHeader() {
  double shift = 0;
  // n_term strict; c_term nostrict; prot n_term match; prot c_term no_match
  // pep n_term no_match; pep c_term no_match
  return DiagonalHeaderPtr(
      new DiagonalHeader(shift, true, false, true, false, false, false));
}

DiagonalHeaderPtr getBottomRightCornerHeader(ProteoformPtr proteoform,
                                             PrmMsPtr ms_six) {
  double prec_mass = ms_six->getHeaderPtr()->getPrecMonoMass();
  double mole_mass = proteoform->getResSeqPtr()->getSeqMass();
  double shift = prec_mass - mole_mass;
  // n term nostrict, c_term strict, prot n_term no match ; prot c_term match
  // pep n_term no match, pep c_term no match 
  return DiagonalHeaderPtr(
      new DiagonalHeader(shift, false, true, false, true, false, false));
}

int findSimilarShiftPos(const std::vector<double> &shifts, double s) {
  int best_pos = -1;
  double best_diff = std::numeric_limits<double>::infinity();
  for(unsigned int i = 0; i < shifts.size();i++){
    if(std::abs(shifts[i] - s) < best_diff){
      best_pos = i;
      best_diff = std::abs(shifts[i] - s);
    }
  }
  return best_pos;
}

bool isExist(const DiagonalHeaderPtrVec &headers, double shift) {
  for(unsigned int i = 0; i < headers.size();i++){
    if(headers[i]->getProtNTermShift() == shift) {
      return true;
    }
  }
  return false;
}

DiagonalHeaderPtrVec PtmSlowMatch::getNTermShiftHeaders(
    std::vector<double> best_shifts,
    PrmMsPtr ms_six,
    ProteoformPtr proteoform,
    PtmMngPtr mng){

  DiagonalHeaderPtrVec headers = getNTermShiftListCommonHeaders(best_shifts);

  DiagonalHeaderPtrVec n_extend_headers;
  // get top-left corner header in spectral grid (shift is 0)
  DiagonalHeaderPtr top_left_corner_header = getTopLeftCornerHeader();
  n_extend_headers.push_back(top_left_corner_header);

  DiagonalHeaderPtrVec c_extend_headers;
  // get bottom-right corner header in the spectral grid. 
  DiagonalHeaderPtr bottom_right_corner_header 
      = getBottomRightCornerHeader(proteoform, ms_six);
  c_extend_headers.push_back(bottom_right_corner_header);

  double prec_mass = ms_six->getHeaderPtr()->getPrecMonoMass();
  std::vector<double> prm_masses = proteoform->getBpSpecPtr()->getPrmMasses();
  // shifts for n_term matches
  std::vector<double> n_term_match_shifts;
  // shifts for c_term matches
  std::vector<double> c_term_match_shifts;
  for(unsigned int i=1;i< prm_masses.size();i++){
    n_term_match_shifts.push_back(-prm_masses[i]);
    c_term_match_shifts.push_back(prec_mass - prm_masses[i]);
  }

  // add trunc headers that have similar shift to best shift headers
  for (unsigned int i = 0; i < headers.size(); i++) {
    double s = headers[i]->getProtNTermShift();
    // find a similar shift in n_term_match_shifts
    int best_n_pos = findSimilarShiftPos(n_term_match_shifts, s);
    if (best_n_pos >= 0) {
      double new_shift = n_term_match_shifts[best_n_pos];
      if (!isExist(n_extend_headers, new_shift)) {
        // n_term strict; c_term nostrict; prot n_term no_match; prot c_term no_match
        // pep n_term match; pep c_term no_match
        DiagonalHeaderPtr header_ptr(
            new DiagonalHeader(new_shift, true, false, false, false, true, false));
        n_extend_headers.push_back(header_ptr);
      }
    }

    int best_c_pos = findSimilarShiftPos(c_term_match_shifts, s);
    if (best_c_pos >= 0) {
      double new_shift = c_term_match_shifts[best_c_pos];
      if (!isExist(c_extend_headers, new_shift)) {
        // n term nostrict, c_term strict, prot n_term no match ; prot c_term no match
        // pep n_term no match, pep c_term match 
        DiagonalHeaderPtr header_ptr(
            new DiagonalHeader(new_shift, false, true, false, false, false, true));
        c_extend_headers.push_back(header_ptr);
      }
    }

  }
  headers.insert(headers.end(), n_extend_headers.begin(), n_extend_headers.end());
  headers.insert(headers.end(), c_extend_headers.begin(), c_extend_headers.end());

  return headers;
}


void PtmSlowMatch::compute(SemiAlignTypePtr type, PrsmPtrVec &prsms) {
  ps_align_->compute(type);
  for (int s = 1; s <= mng_->n_unknown_shift_; s++) {
    PrsmPtr prsm_ptr = geneResult(s);
    prsms.push_back(prsm_ptr);
  }
}


PrsmPtr PtmSlowMatch::geneResult(int shift_num){
  DiagonalHeaderPtrVec headers= ps_align_->getResult(shift_num);
  if(headers.size()==0) {
    return nullptr;
  }
  double refine_prec_mass = ms_three_->getHeaderPtr()->getPrecMonoMass();
  int first_pos = headers[0]->getTruncFirstResPos();
  int last_pos = headers[headers.size()-1]->getTruncLastResPos();
  ProteoformPtr proteoform  = getSubProteoform(proteoform_, first_pos, last_pos);
  DiagonalHeaderPtrVec refined_headers = refineHeadersBgnEnd(
      first_pos, last_pos, proteoform_, deconv_ms_, ms_three_,mng_, headers);

  if(refined_headers.size()==0){
    return nullptr;
  }


  ChangePtrVec changes = getUnexpectedChanges(refined_headers, 
                                              first_pos, last_pos);
  proteoform->addUnexpectedChangePtrVec(changes);

  return PrsmPtr(new Prsm(proteoform, deconv_ms_, refine_prec_mass,
          0, mng_->prsm_para_ptr_->getSpParaPtr()));
}

} /* namespace prot */
