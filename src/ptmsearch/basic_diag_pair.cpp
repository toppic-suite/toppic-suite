/*
 * basic_diag_pair.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

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

BasicDiagPair::BasicDiagPair(int x,int y,double score,
                             int diag_order,double diff,
                             int prm_peak_type,
                             BasicDiagPairDiagPtr diagonal):Pair(x,y){
    score_=0.0;
    diag_order_= diag_order;
    diff_=diff;
    base_type_ = prm_peak_type;
    diagonal_=diagonal;
}

BasicDiagPair::BasicDiagPair(BasicDiagPairPtr pair)
  :Pair(pair->getX(),pair->getY()){
    score_=0.0;
    //can be null in java;
    base_type_=-1;
    diag_order_= pair->getDiagOrder();
    diff_=pair->getDiff();
}

BasicDiagPairPtrVec compDiagPair(PrmMsPtr sp,std::vector<double> seq_masses,
                                 DiagonalHeaderPtr header){
    unsigned int i=1;
    unsigned int j=0;
    std::vector<std::vector<double>> scores;
    std::vector<std::vector<int>> positions;
    for(unsigned int k =0;k< seq_masses.size();k++){
        std::vector<double> score {0,0};
        std::vector<int> position {0,0};
        scores.push_back(score);
        positions.push_back(position);
    }
    double n_term_shift = header->getProtNTermShift();
    std::vector<double> real_masses;
    for(unsigned int k =0;k<sp->size();k++){
        real_masses.push_back(sp->getPeakPtr(k)->getPosition());
    }
    while(i<sp->size()-1 && j<seq_masses.size()){
        PrmPeakPtr peak = sp->getPeakPtr(i);
        int type = peak->getBaseType();
        double error = 0;
        if(header->isNStrict() && !header->isCStrict()){
            error = peak->getNStrictCRelaxTolerance();
        }
        else{
            error = peak->getNRelaxCStrictTolerance();
        }
        double deviation = peak->getPosition()-seq_masses[j] - n_term_shift;
        if(std::abs(deviation) <= error){
            if(scores[j][type]<peak->getScr()){
                scores[j][type] = peak->getScr();
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
    for(j=0;j<seq_masses.size();j++){
        for(int k=0;k<2;k++){
            int pos = positions[j][k];
      if(pos >0){
        double mass = sp->getPeakPtr(pos)->getPosition() - seq_masses[j];
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

    std::sort(pair_list.begin(),pair_list.end(),prot::comparePairUp);
    return pair_list;
}

BasicDiagPairDiagPtrVec getDiagonals(DiagonalHeaderPtrVec headers,
                                     PrmMsPtr ms_six,ProteoformPtr seq,
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
                                 ProteoformPtr seq,
                                 PtmMngPtr mng){
  BpSpecPtr bp_spec_ptr = seq->getBpSpecPtr();
  double n_shift = header->getProtNTermShift();
  double c_shift = ms_six->getHeaderPtr()->getPrecMonoMass()
      -seq->getResSeqPtr()->getSeqMass()-n_shift;
  double term_error_tolerance = mng->term_match_error_tolerance;
  setPrefixSuffix(header,c_shift,seq, term_error_tolerance, mng);
  IonTypePtr b_ion = IonTypeFactory::getIonTypePtr_B();
  std::vector<double> b_masses = bp_spec_ptr->getBreakPointMasses(b_ion);
  BasicDiagPairPtrVec diag_pair_list = compDiagPair(ms_six, b_masses, header);
  if (diag_pair_list.size() > 0) {
    header->setId(cnt);
    BasicDiagPairDiagPtr temp 
        = BasicDiagPairDiagPtr(
            new Diagonal<BasicDiagPairPtr>(header,diag_pair_list));
    for(unsigned int i=0;i<diag_pair_list.size();i++){
      diag_pair_list[i]->setDiagonal(temp);
    }
    return temp;
  }
  return nullptr;
}

bool contains(BasicDiagPairPtrVec pairs,int y){
    for(unsigned int i=0;i<pairs.size();i++){
        if(y==pairs[i]->getY()){
            return true;
        }
    }
    return false;
}
} /* namespace prot */
