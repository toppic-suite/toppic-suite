/*
 * ptm_slow_match.cpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#include <ptmsearch/ptm_slow_match.hpp>
#include <iostream>

namespace prot {

PtmSlowMatch::PtmSlowMatch(ProteoformPtr seq,SpectrumSetPtr spectrum_set,
                           CompShiftLowMemPtr comp_shift,PtmMngPtr mng){
  mng_=mng;
  deconv_ms_ = spectrum_set->getDeconvMs();
  ms_six_=spectrum_set->getSpSix();
  ms_three_ = spectrum_set->getSpThree();
  seq_=seq;
  init(comp_shift);
}

void PtmSlowMatch::init(CompShiftLowMemPtr comp_shift){
  double scale = mng_->ptm_fast_filter_scale_;
  std::vector<std::vector<int>> sp_masses = getIntMassErrorList(ms_six_,scale,true,false);

  std::vector<double> best_shift= comp_shift->findBestShift(
      sp_masses[0],
      sp_masses[1],
      seq_->getBpSpecPtr()->getScaledPrmMasses(scale),
      mng_->n_top_diagonals_,
      mng_->min_diagonal_gap_,
      scale);
  DiagonalHeaderPtrVec n_term_shifts =getNTermShiftList(best_shift,ms_six_,seq_,mng_);
  BasicDiagPairDiagPtrVec diagonals = getDiagonals(n_term_shifts,ms_six_,seq_,mng_);
  std::vector<double> ms_masses = prot::getMassList(ms_six_);
  std::vector<double> seq_masses = seq_->getBpSpecPtr()->getPrmMasses();
  ps_align_ = PSAlignPtr(new PSAlign(ms_masses,seq_masses,diagonals,mng_));
}

void PtmSlowMatch::compute(SemiAlignTypePtr type, PrsmPtrVec &prsms) {
  ps_align_->compute(type);
  scores_ = ps_align_->getAlignScr();
  headers_ = ps_align_->getResult();
  for (int s = 1; s <= mng_->n_unknown_shift_; s++) {
    PrsmPtr prsm_ptr = geneResult(s);
    prsms.push_back(prsm_ptr);
  }
}


PrsmPtr PtmSlowMatch::geneResult(int shift_num){
  DiagonalHeaderPtrVec headers=headers_[shift_num];
  if(headers.size()==0) {
    return nullptr;
  }
  //double refine_prec_mass = ms_three_->getHeaderPtr()->getPrecMonoMass()
  //    +result_deltas_[shift_num];
  //    result_deltas are not used any more 
  double refine_prec_mass = ms_three_->getHeaderPtr()->getPrecMonoMass();
  int first_pos = headers[0]->getTruncFirstResPos();
  int last_pos = headers[headers.size()-1]->getTruncLastResPos();
  DiagonalHeaderPtrVec refined_headers = refineHeadersBgnEnd(
      first_pos, seq_, deconv_ms_, ms_three_,mng_, headers);

  if(refined_headers.size()==0){
    return nullptr;
  }

  ProteoformPtr proteoform  = getSubProteoform(seq_, first_pos, last_pos);

//  std::cout<<refined_headers.size()<<std::endl;
  ChangePtrVec changes = getUnexpectedChanges(
      refined_headers, first_pos, last_pos);
  proteoform->addUnexpectedChangePtrVec(changes);

  return PrsmPtr(
      new Prsm(proteoform, deconv_ms_, refine_prec_mass,
          0, mng_->prsm_para_ptr_->getSpParaPtr()));
}

DiagonalHeaderPtrVec PtmSlowMatch::getNTermShiftList(
    std::vector<double> best_shift,
    PrmMsPtr ms_six,
    ProteoformPtr seq,
    PtmMngPtr mng){

  DiagonalHeaderPtrVec headers = getNTermShiftListCommon(best_shift);
  // get protein N term shift 
  DiagonalHeaderPtrVec n_term_shifts_comp_left = getNTermShiftListCompLeft(seq,mng);
  // get protein C term shift 
  DiagonalHeaderPtrVec n_term_shifts_comp_right = getNTermShiftListCompRight(seq,ms_six);
  DiagonalHeaderPtrVec extend_n_term_shifts;

  for(unsigned int i=0;i<n_term_shifts_comp_left.size();i++){
    headers.push_back(n_term_shifts_comp_left[i]);
  }
  for(unsigned int i=0;i<n_term_shifts_comp_right.size();i++){
    headers.push_back(n_term_shifts_comp_right[i]);
  }
  std::vector<double> ms_masses = getMassList(ms_six);
  std::vector<double> seq_masses = seq->getBpSpecPtr()->getPrmMasses();
  double shift;
  // if a diagonal has a similar shift to a trunc diagonal, then trunc diagonal
  // is added to the list.
  for(unsigned int i=1;i<seq_masses.size();i++){
    shift = - seq_masses[i];
    if(found(shift,headers,mng_->extend_trunc_error_tolerance_)){
      extend_n_term_shifts.push_back(DiagonalHeaderPtr(
              new DiagonalHeader(shift,true,false,true,false)));
    }
  }

  for(unsigned int i=1;i<seq_masses.size();i++){
    shift = ms_masses[ms_masses.size()-1] - seq_masses[i];
    if(found(shift,headers,mng_->extend_trunc_error_tolerance_)){
      extend_n_term_shifts.push_back(DiagonalHeaderPtr(
              new DiagonalHeader(shift,false,true,false,true)));
    }
  }

  for(unsigned int i=0;i<extend_n_term_shifts.size();i++){
    headers.push_back(extend_n_term_shifts[i]);
  }

  return headers;
}

bool PtmSlowMatch::found(double shift, DiagonalHeaderPtrVec headerlist, 
                         double trunc_error_tolerance){
  for(unsigned int i=0;i<headerlist.size();i++){
    if(std::abs(shift-headerlist[i]->getProtNTermShift())
       <= trunc_error_tolerance){
      return true;
    }
  }
  return false;
}
} /* namespace prot */
