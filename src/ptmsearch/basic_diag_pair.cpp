#include "ptmsearch/basic_diag_pair.hpp"

namespace prot {

BasicDiagPair::BasicDiagPair(int x,int y,double score,
                             int diag_order,double diff)
    :Pair(x,y){

    score_ = score;
    diag_order_= diag_order;
    diff_=diff;
}

inline BasicDiagPairPtrVec compDiagPair(const PrmPeakPtrVec &prm_peaks, 
                                        int group_spec_num,  
                                        const std::vector<double> &seq_masses,
                                        DiagonalHeaderPtr header_ptr){
  std::vector<std::vector<double>> n_term_scores(
      group_spec_num, std::vector<double>(seq_masses.size(), 0));
  std::vector<std::vector<double>> c_term_scores(
      group_spec_num, std::vector<double>(seq_masses.size(), 0));
  std::vector<double> top_scores(seq_masses.size(), 0);
  std::vector<int> top_positions(seq_masses.size(), -1);

  // i starts from 0 and ends at size - 1 to include the first prm 0 and the 
  // last prm precursor_mass - water_mass
  size_t i=0;
  size_t j=0;
  double n_term_shift = header_ptr->getProtNTermShift();
  std::vector<double> real_masses;
  for(size_t k =0;k< prm_peaks.size();k++){
    real_masses.push_back(prm_peaks[k]->getPosition());
  }
  //LOG_DEBUG("real mass size " << real_masses.size() << " seq mass size " << seq_masses.size());
  while(i < prm_peaks.size() && j < seq_masses.size()){
    PrmPeakPtr peak = prm_peaks[i];
    int type = peak->getBaseType();
    int spec_id = peak->getSpecId();
    double error = 0;
    if(header_ptr->isNStrict() && !header_ptr->isCStrict()){
      error = peak->getNStrictCRelaxTolerance();
    }
    else{
      error = peak->getNRelaxCStrictTolerance();
    }
    double deviation = peak->getPosition()-seq_masses[j] - n_term_shift;
    //LOG_DEBUG("deviation" << i << " " << j << " spec id " << spec_id);
    if(std::abs(deviation) <= error){
      double peak_score = peak->getScore();
      if (type == PRM_PEAK_TYPE_ORIGINAL) {
        if(n_term_scores[spec_id][j] < peak_score){
          n_term_scores[spec_id][j] = peak_score;
        }
      }
      else {
        if(c_term_scores[spec_id][j] < peak_score){
          c_term_scores[spec_id][j] = peak_score;
        }
      }
      // update top position
      if (top_scores[j] < peak_score) {
        top_scores[j] = peak_score;
        top_positions[j] = i;
      }
    }
    //LOG_DEBUG("start increase  peak " << peak.get());
    if(increaseIJ(i,j,deviation,peak->getNRelaxCStrictTolerance(),
                  real_masses,seq_masses)){
      i++;
    }
    else{
      j++;
    }
    //LOG_DEBUG("increase" << i << " " << j);
  }
  std::vector<double> sum_scores(seq_masses.size(), 0);
  for (size_t p = 0; p < seq_masses.size(); p++) {
    for (int m = 0; m < group_spec_num; m++) {
      sum_scores[p] += n_term_scores[m][p];
      sum_scores[p] += c_term_scores[m][p];
    }
  }

  // add pairs 
  BasicDiagPairPtrVec  pair_list;
  for(j=0; j<seq_masses.size();j++){
    int pos = top_positions[j];
    if(pos >= 0){
      double diff = prm_peaks[pos]->getPosition() - seq_masses[j];
      double score = sum_scores[j];
      BasicDiagPairPtr diag_pair_ptr(
          new BasicDiagPair(pos,j, score ,pair_list.size(),diff));
      pair_list.push_back(diag_pair_ptr);
    }
  }
  // compare y (seq_mass) first, then x (peak position)
  std::sort(pair_list.begin(),pair_list.end(),comparePairUp);
  return pair_list;
}

BasicDiagonalPtr getDiagonalPtr(DiagonalHeaderPtr header_ptr,
                                const PrmPeakPtrVec &prm_peaks,
                                double prec_mono_mass,
                                int group_spec_num,
                                ProteoformPtr proteo_ptr,
                                PtmMngPtr mng_ptr){
  BpSpecPtr bp_spec_ptr = proteo_ptr->getBpSpecPtr();
  double n_shift = header_ptr->getProtNTermShift();
  double c_shift = prec_mono_mass - proteo_ptr->getResSeqPtr()->getSeqMass() - n_shift;
  header_ptr->initData(c_shift, proteo_ptr, mng_ptr->align_prefix_suffix_shift_thresh_);

  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  BasicDiagPairPtrVec diag_pair_list = compDiagPair(prm_peaks, group_spec_num, prm_masses, header_ptr);
  if (diag_pair_list.size() > 0) {
    BasicDiagonalPtr diagonal_ptr 
        = BasicDiagonalPtr(new Diagonal<BasicDiagPairPtr>(header_ptr,diag_pair_list));
    for(size_t i=0;i<diag_pair_list.size();i++){
      diag_pair_list[i]->setDiagonalPtr(diagonal_ptr);
    }
    return diagonal_ptr;
  }
  return nullptr;
}

BasicDiagonalPtrVec getDiagonals(const DiagonalHeaderPtrVec& header_ptr_vec,
                                 const PrmPeakPtrVec &prm_peaks, 
                                 double prec_mono_mass, int group_spec_num,
                                 ProteoformPtr proteo_ptr, PtmMngPtr mng_ptr){
  BasicDiagonalPtrVec diagonal_list;
  for(size_t i=0; i<header_ptr_vec.size(); i++){
    BasicDiagonalPtr diagonal_ptr 
        = getDiagonalPtr(header_ptr_vec[i], prm_peaks, prec_mono_mass, 
                         group_spec_num, proteo_ptr, mng_ptr);
    if(diagonal_ptr!=nullptr){
      diagonal_list.push_back(diagonal_ptr);
    }
  }
  // important set id for headers
  for (size_t i = 0; i < diagonal_list.size(); i++) {
    diagonal_list[i]->getHeader()->setId(i);
  }
  return diagonal_list;
}


/*
inline BasicDiagPairPtrVec compDiagPair(PrmMsPtr ms_ptr, 
                                        const std::vector<double>& seq_masses,
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
  // add pairs 
  BasicDiagPairPtrVec  pair_list;
  for(j=0; j<seq_masses.size();j++){
    int pos = n_term_positions[j];
    if(pos >= 0){
      double diff = ms_ptr->getPeakPtr(pos)->getPosition() - seq_masses[j];
      BasicDiagPairPtr diag_pair_ptr(new BasicDiagPair(pos,j,1.0,pair_list.size(),
                                                       diff, PRM_PEAK_TYPE_ORIGINAL));
      pair_list.push_back(diag_pair_ptr);
    }
    pos = c_term_positions[j];
    if (pos >= 0) {
      double diff = ms_ptr->getPeakPtr(pos)->getPosition() - seq_masses[j];
      BasicDiagPairPtr diag_pair_ptr(new BasicDiagPair(pos,j,1.0,pair_list.size(),
                                                       diff, PRM_PEAK_TYPE_REVERSED));
      pair_list.push_back(diag_pair_ptr);
    }
  }

  // compare y (seq_mass) first, then x (peak position)
  std::sort(pair_list.begin(),pair_list.end(),comparePairUp);
  return pair_list;
}

BasicDiagPairPtrVec compDiagPair(const PrmMsPtrVec &ms_ptr_vec, 
                                 const std::vector<double>& seq_masses,
                                 DiagonalHeaderPtr header_ptr) {
  BasicDiagPairPtrVec  pair_list;
  for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
    BasicDiagPairPtrVec pairs = compDiagPair(ms_ptr_vec[i], seq_masses, header_ptr);
    pair_list.insert(pair_list.end(), pairs.begin(), pairs.end());
  }
  return pair_list;
}
*/



} /* namespace prot */
