// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "oneptmsearch/basic_diag_pair.hpp"

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
    BasePeakTypePtr type_ptr = peak->getBaseTypePtr();
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
      if (type_ptr == BasePeakType::ORIGINAL) {
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
  std::sort(pair_list.begin(),pair_list.end(),Pair::cmpPosInc);
  return pair_list;
}

BasicDiagonalPtr getDiagonalPtr(DiagonalHeaderPtr header_ptr,
                                const PrmPeakPtrVec &prm_peaks,
                                int group_spec_num,
                                ProteoformPtr proteo_ptr){
  BpSpecPtr bp_spec_ptr = proteo_ptr->getBpSpecPtr();
  /*double n_shift = header_ptr->getProtNTermShift();*/
  //double c_shift = prec_mono_mass - proteo_ptr->getResSeqPtr()->getSeqMass() - n_shift;
  /*header_ptr->initData(c_shift, proteo_ptr, mng_ptr->align_prefix_suffix_shift_thresh_);*/

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

BasicDiagonalPtrVec geneDiagonals(const DiagonalHeaderPtrVec& header_ptr_vec,
        const PrmPeakPtrVec &prm_peaks, 
        int group_spec_num, ProteoformPtr proteo_ptr){

    BasicDiagonalPtrVec diagonal_list;
    for(size_t i=0; i<header_ptr_vec.size(); i++){
        BasicDiagonalPtr diagonal_ptr 
            = getDiagonalPtr(header_ptr_vec[i], prm_peaks, 
                    group_spec_num, proteo_ptr);
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

} /* namespace prot */
