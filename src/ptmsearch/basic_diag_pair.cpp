#include "ptmsearch/basic_diag_pair.hpp"

namespace prot {

BasicDiagPair::BasicDiagPair(int x,int y,double score,
                             int diag_order,double diff,
                             int prm_peak_type):Pair(x,y){
    score_ = score;
    diag_order_= diag_order;
    diff_=diff;
    base_type_ = prm_peak_type;
}

BasicDiagPairPtrVec compDiagPair(PrmMsPtr ms_ptr, const std::vector<double>& seq_masses,
                                 DiagonalHeaderPtr header_ptr){
  std::vector<double> n_term_scores(seq_masses.size());
  std::vector<double> c_term_scores(seq_masses.size());
  std::vector<int> n_term_positions(seq_masses.size(), -1);
  std::vector<int> c_term_positions(seq_masses.size(), -1);

  // i starts from 0 and ends at size - 1 to include the first prm 0 and the 
  // last prm precursor_mass - water_mass
  size_t i=0;
  size_t j=0;
  double n_term_shift = header_ptr->getProtNTermShift();
  std::vector<double> real_masses;
  for(size_t k =0;k<ms_ptr->size();k++){
    real_masses.push_back(ms_ptr->getPeakPtr(k)->getPosition());
  }

  while(i < ms_ptr->size() && j < seq_masses.size()){
    PrmPeakPtr peak = ms_ptr->getPeakPtr(i);
    int type = peak->getBaseType();
    double error = 0;
    if(header_ptr->isNStrict() && !header_ptr->isCStrict()){
      error = peak->getNStrictCRelaxTolerance();
    }
    else{
      error = peak->getNRelaxCStrictTolerance();
    }
    double deviation = peak->getPosition()-seq_masses[j] - n_term_shift;
    if(std::abs(deviation) <= error){
      if (type == PRM_PEAK_TYPE_ORIGINAL) {
        if(n_term_scores[j] < peak->getScore()){
          n_term_scores[j] = peak->getScore();
          n_term_positions[j] = i;
        }
      }
      else {
        if(c_term_scores[j] < peak->getScore()){
          c_term_scores[j] = peak->getScore();
          c_term_positions[j] = i;
        }
      }
    }
    if(increaseIJ(i,j,deviation,peak->getNRelaxCStrictTolerance(),
                  real_masses,seq_masses)){
      i++;
    }
    else{
      j++;
    }
  }
  /* add pairs */
  BasicDiagPairPtrVec pair_list;
  for(j=0; j<seq_masses.size();j++){
    int pos = n_term_positions[j];
    if(pos >= 0){
      double diff = ms_ptr->getPeakPtr(pos)->getPosition() - seq_masses[j];
      pair_list.push_back(BasicDiagPairPtr(
              new BasicDiagPair(pos,j,1.0,pair_list.size(),
                                diff, PRM_PEAK_TYPE_ORIGINAL)));
    }
    pos = c_term_positions[j];
    if (pos >= 0) {
      double diff = ms_ptr->getPeakPtr(pos)->getPosition() - seq_masses[j];
      pair_list.push_back(BasicDiagPairPtr(
              new BasicDiagPair(pos,j,1.0,pair_list.size(),
                                diff, PRM_PEAK_TYPE_REVERSED)));
    }
  }

  // compare y (seq_mass) first, then x (peak position)
  std::sort(pair_list.begin(),pair_list.end(),comparePairUp);
  return pair_list;
}

BasicDiagPairDiagPtrVec getDiagonals(const DiagonalHeaderPtrVec& header_ptrs,
                                     PrmMsPtr ms_six_ptr, ProteoformPtr proteo_ptr,
                                     PtmMngPtr mng_ptr){
  BasicDiagPairDiagPtrVec diagonal_list;
  int cnt =0;
  for(size_t i=0;i<header_ptrs.size();i++){
    BasicDiagPairDiagPtr diagonal_ptr = getDiagonal(cnt,header_ptrs[i],ms_six_ptr,
                                                    proteo_ptr, mng_ptr);
    if(diagonal_ptr!=nullptr){
      diagonal_list.push_back(diagonal_ptr);
      cnt++;
    }
  }
  return diagonal_list;
}

BasicDiagPairDiagPtr getDiagonal(int cnt,DiagonalHeaderPtr header_ptr,
                                 PrmMsPtr ms_six_ptr,
                                 ProteoformPtr proteo_ptr,
                                 PtmMngPtr mng_ptr){
  BpSpecPtr bp_spec_ptr = proteo_ptr->getBpSpecPtr();
  double n_shift = header_ptr->getProtNTermShift();
  double c_shift = ms_six_ptr->getHeaderPtr()->getPrecMonoMass()
      -proteo_ptr->getResSeqPtr()->getSeqMass()-n_shift;
  header_ptr->initData(c_shift, proteo_ptr, mng_ptr->align_prefix_suffix_shift_thresh_);

  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  BasicDiagPairPtrVec diag_pair_list = compDiagPair(ms_six_ptr, prm_masses, header_ptr);
  if (diag_pair_list.size() > 0) {
    header_ptr->setId(cnt);
    BasicDiagPairDiagPtr diagonal_ptr 
        = BasicDiagPairDiagPtr(new Diagonal<BasicDiagPairPtr>(header_ptr,diag_pair_list));
    for(size_t i=0;i<diag_pair_list.size();i++){
      diag_pair_list[i]->setDiagonal(diagonal_ptr);
    }
    return diagonal_ptr;
  }
  return nullptr;
}

} /* namespace prot */
