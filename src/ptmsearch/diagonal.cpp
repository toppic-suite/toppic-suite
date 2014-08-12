
#include "base/logger.hpp"
#include "ptmsearch/diagonal.hpp"

namespace prot {

inline TheoPeakPtrVec getDiagonalTheoPeak(ProteoformPtr proteo_ptr, 
                                          ActivationPtr activation_ptr,
                                          const DiagonalHeaderPtrVec &header_ptrs,
                                          int i, double min_mass){
  DiagonalHeaderPtr first_header_ptr = header_ptrs[0];
  DiagonalHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size()-1];
  int first_res_pos = first_header_ptr->getTruncFirstResPos();
  int last_res_pos = last_header_ptr->getTruncLastResPos();
  ResSeqPtr subseq_ptr = proteo_ptr->getResSeqPtr()->getSubResidueSeq(first_res_pos,last_res_pos);
  BpSpecPtr pep_ptr = BpSpecPtr(new BpSpec(subseq_ptr));
  double pep_n_term_shift = header_ptrs[i]->getProtNTermShift()
      -first_header_ptr->getProtNTermShift()
      +first_header_ptr->getPepNTermShift();
  double pep_c_term_shift = header_ptrs[i]->getProtCTermShift()
      -last_header_ptr->getProtCTermShift()
      +last_header_ptr->getPepCTermShift();
  double max_mass = subseq_ptr->getSeqMass() + 
      pep_n_term_shift + pep_c_term_shift - min_mass;
  return getTheoPeak(pep_ptr,activation_ptr,NeutralLossFactory::getNeutralLossPtr_NONE(),
                     pep_n_term_shift,pep_c_term_shift,
                     header_ptrs[i]->getMatchFirstBpPos()-first_res_pos,
                     header_ptrs[i]->getMatchLastBpPos()-first_res_pos,
                     min_mass, 
                     max_mass);
}


DiagonalHeaderPtrVec refineHeadersBgnEnd(
    int first_res_pos,
    int last_res_pos,
    ProteoformPtr proteo_ptr,
    DeconvMsPtr deconv_ms_ptr,
    ExtendMsPtr ms_three_ptr,
    PtmMngPtr mng_ptr,
    const DiagonalHeaderPtrVec& header_ptrs){

  DiagonalHeaderPtrVec result_list;
  for(size_t i=0;i<header_ptrs.size();i++){
    TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
        proteo_ptr,
        deconv_ms_ptr->getHeaderPtr()->getActivationPtr(),
        header_ptrs,
        i,
        mng_ptr->prsm_para_ptr_->getSpParaPtr()->getMinMass());

    int bgn = header_ptrs[i]->getMatchFirstBpPos()-first_res_pos;
    int end = header_ptrs[i]->getMatchLastBpPos()-first_res_pos;
    PeakIonPairPtrVec pair_ptrs = findPairs(ms_three_ptr, theo_peak_ptrs, bgn, end);
    if(pair_ptrs.size()<1){
      int pair_size = pair_ptrs.size();
      LOG_WARN("Empty Segment is found "+prot::convertToString(pair_size));
      if (i == 0 ) {
        int new_bgn = first_res_pos;
        int new_end = first_res_pos;
        header_ptrs[i]->setMatchFirstBpPos(new_bgn);
        header_ptrs[i]->setMatchLastBpPos(new_end);
        result_list.push_back(header_ptrs[i]);
      }
      else if (i == header_ptrs.size() - 1) {
        int new_bgn = last_res_pos + 1;
        int new_end = last_res_pos + 1;
        header_ptrs[i]->setMatchFirstBpPos(new_bgn);
        header_ptrs[i]->setMatchLastBpPos(new_end);
        result_list.push_back(header_ptrs[i]);
      }
    }
    else{
      int new_bgn = first_res_pos + getNewBgn(pair_ptrs);
      int new_end = first_res_pos + getNewEnd(pair_ptrs);
      header_ptrs[i]->setMatchFirstBpPos(new_bgn);
      header_ptrs[i]->setMatchLastBpPos(new_end);
      result_list.push_back(header_ptrs[i]);
    }
  }
  return result_list;
}

int getNewBgn(const PeakIonPairPtrVec &pair_ptrs){
  int new_bgn = std::numeric_limits<int>::max();
  for(size_t i=0;i<pair_ptrs.size();i++) {
    if(pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos() < new_bgn){
      new_bgn = pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
    }
  }
  return new_bgn;
}
int getNewEnd(const PeakIonPairPtrVec &pair_ptrs){
  int new_end = 0;
  for(size_t i=0;i<pair_ptrs.size();i++) {
    if(pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos() > new_end){
      new_end = pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
    }
  }
  return new_end;
}

} /* namespace prot */
