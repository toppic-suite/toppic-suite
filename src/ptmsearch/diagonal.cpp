/*
 * diagonal.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#include "ptmsearch/diagonal.hpp"
#include "base/logger.hpp"
#include "float.h"
#include "limits.h"

namespace prot {

TheoPeakPtrVec getDiagonalTheoPeak(ProteoformPtr seq,ActivationPtr type,
                                   DiagonalHeaderPtrVec headers,int i,
                                   double minMass){
  DiagonalHeaderPtr first_header = headers[0];
  DiagonalHeaderPtr last_header = headers[headers.size()-1];
  int first_res_pos = first_header->getTruncFirstResPos();
  int last_res_pos = last_header->getTruncLastResPos();
  ResSeqPtr subseq = seq->getResSeqPtr()->getSubResidueSeq(first_res_pos,last_res_pos);
  BpSpecPtr pep = BpSpecPtr(new BpSpec(subseq));
  double pep_n_term_shift = headers[i]->getProtNTermShift()
      -first_header->getProtNTermShift()
      +first_header->getPepNTermShift();
  double pep_c_term_shift = headers[i]->getProtCTermShift()
      -last_header->getProtCTermShift()
      +last_header->getPepCTermShift();
  double maxMass = subseq->getSeqMass() + 
      pep_n_term_shift + pep_c_term_shift - minMass;
  //  std::cout<<std::fixed<<pep->getResSeqMass()<<"|"<<pep_n_term_shift<<"|"<<pep_c_term_shift<<std::endl;
  return getTheoPeak(pep,type,NeutralLossFactory::getNeutralLossPtr_NONE(),
                     pep_n_term_shift,pep_c_term_shift,
                     headers[i]->getMatchFirstBpPos()-first_res_pos,
                     headers[i]->getMatchLastBpPos()-first_res_pos,
                     minMass, 
                     maxMass);
}

DiagonalHeaderPtrVec refineHeadersBgnEnd(
    int first_res_pos,
    int last_res_pos,
    ProteoformPtr seq,
    DeconvMsPtr deconv_ms,
    ExtendMsPtr ms_three,
    PtmMngPtr mng,
    DiagonalHeaderPtrVec headers){
  DiagonalHeaderPtrVec result_list;
  //    std::cout<<headers.size()<<std::endl;
  for(unsigned int i=0;i<headers.size();i++){
    TheoPeakPtrVec ions = prot::getDiagonalTheoPeak(
        seq,
        deconv_ms->getHeaderPtr()->getActivationPtr(),
        headers,
        i,
        mng->prsm_para_ptr_->getSpParaPtr()->getMinMass());
    int bgn = headers[i]->getMatchFirstBpPos()-first_res_pos;
    int end = headers[i]->getMatchLastBpPos()-first_res_pos;
    PeakIonPairPtrVec pairs = findPairs(ms_three,ions,bgn,end);
    if(pairs.size()<1){
      int pair_size = pairs.size();
      LOG_WARN("Empty Segment is found "+prot::convertToString(pair_size));
      if (i == 0 ) {
        int new_bgn = first_res_pos;
        int new_end = first_res_pos;
        headers[i]->setMatchFirstBpPos(new_bgn);
        headers[i]->setMatchLastBpPos(new_end);
        result_list.push_back(headers[i]);
      }
      else if (i == headers.size() - 1) {
        int new_bgn = last_res_pos + 1;
        int new_end = last_res_pos + 1;
        headers[i]->setMatchFirstBpPos(new_bgn);
        headers[i]->setMatchLastBpPos(new_end);
        result_list.push_back(headers[i]);
      }
    }
    else{
      int new_bgn = first_res_pos + getNewBgn(pairs);
      int new_end = first_res_pos + getNewEnd(pairs);
      headers[i]->setMatchFirstBpPos(new_bgn);
      headers[i]->setMatchLastBpPos(new_end);
      result_list.push_back(headers[i]);
    }
  }
  return result_list;
}

int getNewBgn(PeakIonPairPtrVec pairs){
  int new_bgn = std::numeric_limits<int>::max();
  for(unsigned int i=0;i<pairs.size();i++) {
    if(pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos() < new_bgn){
      new_bgn = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
    }
  }
  return new_bgn;
}
int getNewEnd(PeakIonPairPtrVec pairs){
  int new_end = 0;
  for(unsigned int i=0;i<pairs.size();i++) {
    if(pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos() > new_end){
      new_end = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
    }
  }
  return new_end;
}
} /* namespace prot */
