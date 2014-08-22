
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

double refinePrecursorAndHeaderShift(ProteoformPtr proteo_ptr,
                                     ExtendMsPtr ms_three_ptr, 
                                     DiagonalHeaderPtrVec &header_ptrs,
                                     PtmMngPtr mng_ptr) {
  double prec_mass = ms_three_ptr->getHeaderPtr()->getPrecMonoMass();
  if (header_ptrs.size() == 0) {
    return prec_mass;
  }
  double tole = ms_three_ptr->getHeaderPtr()->getErrorTolerance();
  int one_side_step_num = 0;
  if (tole > 0) {
    one_side_step_num = (int)std::floor(tole/mng_ptr->refine_prec_step_width_);
  }
  int step_num = 1 + 2 * one_side_step_num;
  DiagonalHeaderPtrVec test_ptrs;
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    test_ptrs.push_back(header_ptrs[i]->clone());
  }
  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    int middle = (test_ptrs[i]->getMatchLastBpPos() + test_ptrs[i+1]->getMatchFirstBpPos())/2;
    test_ptrs[i]->setMatchLastBpPos(middle);
    test_ptrs[i+1]->setMatchFirstBpPos(middle);
  }
  test_ptrs[0]->setMatchFirstBpPos(test_ptrs[0]->getTruncFirstResPos());
  DiagonalHeaderPtr last_test_ptr = test_ptrs[test_ptrs.size()-1];
  last_test_ptr->setMatchLastBpPos(last_test_ptr->getTruncLastResPos());
  //LOG_DEBUG("TEST PTR SIZE" << test_ptrs.size());
  double change = - one_side_step_num * mng_ptr->refine_prec_step_width_;
  /* for the first test_ptrs.size() - 1 headers, change C-term shift,
   * for the last header, change N-term shift */
  for (size_t i = 0; i < test_ptrs.size() - 1; i++) {
    test_ptrs[i]->changeOnlyCTermShift(change);
  }
  last_test_ptr->changeOnlyNTermShift(change);
  int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
  double best_score = -1;
  int best_i_start = -1;
  int best_i_end = -1;
  bool continuous = false;
  /*
  for (size_t i = 0; i < ms_three_ptr->size(); i++) {
    std::cout << "ms " << ms_three_ptr->getPeakPtr(i)->getPosition() << std::endl;
  }
  */
  for (int i = 0; i < step_num; i++) {
    //double delta = (i - one_side_step_num) * mng_ptr->refine_prec_step_width_;
    double cur_score = 0;
    for(size_t j=0; j< test_ptrs.size(); j++){
      TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
          proteo_ptr,
          ms_three_ptr->getHeaderPtr()->getActivationPtr(),
          test_ptrs,
          j,
          mng_ptr->prsm_para_ptr_->getSpParaPtr()->getMinMass());

      int bgn = test_ptrs[j]->getMatchFirstBpPos()-first_res_pos;
      int end = test_ptrs[j]->getMatchLastBpPos()-first_res_pos;
      PeakIonPairPtrVec pair_ptrs = findPairs(ms_three_ptr, theo_peak_ptrs, bgn, end);
      /*
      std::cout << j << " begin " << bgn << " end " << end << std::endl;
      if (i == 27) {
        for (size_t k = 0; k < theo_peak_ptrs.size(); k++) {
          std::cout << j << " theo " << theo_peak_ptrs[k]->getPosition() << std::endl;
        }
        for (size_t k = 0; k < pair_ptrs.size(); k++) {
          IonPtr ion_ptr = pair_ptrs[k]->getTheoPeakPtr()->getIonPtr();
          std::cout << j << " match " << k << " " << pair_ptrs[k]->getRealPeakPtr()->getPosition() 
              << " " << pair_ptrs[k]->getTheoPeakPtr()->getPosition() << " " 
              << ion_ptr->getIonTypePtr()->getName() << " " 
              << ion_ptr->getDisplayPos() << std::endl;
        }
      }
      */
      cur_score += pair_ptrs.size();
    }
    //LOG_DEBUG("i " << i << " header number "<<  test_ptrs.size() << " delta " << delta << " cur score " << cur_score);
    if (cur_score > best_score) {
      best_score = cur_score;
      best_i_start = i;
      best_i_end = i;
      continuous = true;
    }
    else if (cur_score == best_score) {
      if (continuous) {
        best_i_end = i;
      }
    }
    else {
      continuous = false;
    }

    change = mng_ptr->refine_prec_step_width_;
    for (size_t j = 0; j < test_ptrs.size() - 1; j++) {
      test_ptrs[j]->changeOnlyCTermShift(change);
    }
    last_test_ptr->changeOnlyNTermShift(change);
  }
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    header_ptrs[i]->setMatchFirstBpPos(test_ptrs[i]->getMatchFirstBpPos());
    header_ptrs[i]->setMatchLastBpPos(test_ptrs[i]->getMatchLastBpPos());
  }
  double best_i = (best_i_end + best_i_start) /2; 
  double best_delta = (best_i - one_side_step_num) * mng_ptr->refine_prec_step_width_;
  //LOG_DEBUG("best_delta " << best_delta);
  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    header_ptrs[i]->changeOnlyCTermShift(best_delta);
  }
  DiagonalHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size()-1];
  last_header_ptr->changeOnlyNTermShift(best_delta);

  return prec_mass + best_delta;
}


DiagonalHeaderPtrVec refineHeadersBgnEnd(
    ProteoformPtr proteo_ptr,
    ExtendMsPtr ms_three_ptr,
    const DiagonalHeaderPtrVec& header_ptrs,
    PtmMngPtr mng_ptr){

  DiagonalHeaderPtrVec result_list;
  int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
  int last_res_pos = header_ptrs[header_ptrs.size()-1]->getTruncLastResPos();
  for(size_t i=0;i<header_ptrs.size();i++){
    TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
        proteo_ptr,
        ms_three_ptr->getHeaderPtr()->getActivationPtr(),
        header_ptrs,
        i,
        mng_ptr->prsm_para_ptr_->getSpParaPtr()->getMinMass());

    int bgn = header_ptrs[i]->getMatchFirstBpPos()-first_res_pos;
    int end = header_ptrs[i]->getMatchLastBpPos()-first_res_pos;

    PeakIonPairPtrVec pair_ptrs = findPairs(ms_three_ptr, theo_peak_ptrs, bgn, end);
    if(pair_ptrs.size()<1){
      int pair_size = pair_ptrs.size();
      LOG_TRACE("Empty Segment is found "+prot::convertToString(pair_size));
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
