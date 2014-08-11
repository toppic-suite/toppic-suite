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
  // i starts from 0 and ends at size - 1 to include the first prm 0 and the 
  // last prm precursor_mass - water_mass
  std::vector<std::vector<double>> scores;
  std::vector<std::vector<int>> positions;
  for(unsigned int k =0;k< seq_masses.size();k++){
    std::vector<double> score {0,0};
    std::vector<int> position {-1,-1};
    scores.push_back(score);
    positions.push_back(position);
  }
  unsigned int i=0;
  unsigned int j=0;
  double n_term_shift = header_ptr->getProtNTermShift();
  std::vector<double> real_masses;
  for(unsigned int k =0;k<ms_ptr->size();k++){
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
      if(scores[j][type]<peak->getScore()){
        scores[j][type] = peak->getScore();
        positions[j][type]=i;
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
  BasicDiagPairPtrVec pair_list;
  for(j=0; j<seq_masses.size();j++){
    for(int k=0;k<2;k++){
      int pos = positions[j][k];
      if(pos >= 0){
        double mass = ms_ptr->getPeakPtr(pos)->getPosition() - seq_masses[j];
        if(k==0){
          pair_list.push_back(BasicDiagPairPtr(
                  new BasicDiagPair(pos,j,1.0,pair_list.size(),
                                    mass,
                                    PRM_PEAK_TYPE_ORIGINAL)));
        }
        else{
          pair_list.push_back(BasicDiagPairPtr(
                  new BasicDiagPair(pos,j,1.0,pair_list.size(),
                                    mass, 
                                    PRM_PEAK_TYPE_REVERSED)));
        }
      }
    }
  }

  // compare y (seq_mass) first, then x (peak position)
  std::sort(pair_list.begin(),pair_list.end(),comparePairUp);
  return pair_list;
}

BasicDiagPairDiagPtrVec getDiagonals(const DiagonalHeaderPtrVec& headers,
                                     PrmMsPtr ms_six, ProteoformPtr seq,
                                     PtmMngPtr mng){
  BasicDiagPairDiagPtrVec diagonal_list;
  int cnt =0;
  for(unsigned int i=0;i<headers.size();i++){
    BasicDiagPairDiagPtr diagonal = getDiagonal(cnt,headers[i],ms_six,seq,mng);
    if(diagonal!=nullptr){
      diagonal_list.push_back(diagonal);
      cnt++;
    }
  }
  return diagonal_list;
}

BasicDiagPairDiagPtr getDiagonal(int cnt,DiagonalHeaderPtr header,
                                 PrmMsPtr ms_six,
                                 ProteoformPtr proteoform,
                                 PtmMngPtr mng){
  BpSpecPtr bp_spec_ptr = proteoform->getBpSpecPtr();
  double n_shift = header->getProtNTermShift();
  double c_shift = ms_six->getHeaderPtr()->getPrecMonoMass()
      -proteoform->getResSeqPtr()->getSeqMass()-n_shift;
  header->initData(c_shift, proteoform, mng->align_prefix_suffix_shift_thresh_);

  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  BasicDiagPairPtrVec diag_pair_list = compDiagPair(ms_six, prm_masses, header);
  if (diag_pair_list.size() > 0) {
    header->setId(cnt);
    BasicDiagPairDiagPtr diagonal 
        = BasicDiagPairDiagPtr(new Diagonal<BasicDiagPairPtr>(header,diag_pair_list));
    for(unsigned int i=0;i<diag_pair_list.size();i++){
      diag_pair_list[i]->setDiagonal(diagonal);
    }
    return diagonal;
  }
  return nullptr;
}

} /* namespace prot */
