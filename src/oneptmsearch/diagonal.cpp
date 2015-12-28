#include "base/logger.hpp"
#include "base/neutral_loss_base.hpp"
#include "spec/theo_peak_factory.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "oneptmsearch/diagonal.hpp"

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
  return TheoPeakFactory::geneTheoPeak(pep_ptr,activation_ptr,NeutralLossBase::getNeutralLossPtr_NONE(),
                                       pep_n_term_shift,pep_c_term_shift,
                                       header_ptrs[i]->getMatchFirstBpPos()-first_res_pos,
                                       header_ptrs[i]->getMatchLastBpPos()-first_res_pos,
                                       min_mass, 
                                       max_mass);
}

inline TheoPeakPtrVec getCTermTheoPeakPtrs (const TheoPeakPtrVec theo_peak_ptrs) {
  TheoPeakPtrVec c_term_peak_ptrs; 
  for (size_t i = 0; i < theo_peak_ptrs.size(); i++) {
    if (!theo_peak_ptrs[i]->getIonPtr()->getIonTypePtr()->isNTerm()) {
      c_term_peak_ptrs.push_back(theo_peak_ptrs[i]);
    }
  }
  return c_term_peak_ptrs;
}

inline TheoPeakPtrVec getNTermTheoPeakPtrs (const TheoPeakPtrVec theo_peak_ptrs) {
  TheoPeakPtrVec n_term_peak_ptrs; 
  for (size_t i = 0; i < theo_peak_ptrs.size(); i++) {
    if (theo_peak_ptrs[i]->getIonPtr()->getIonTypePtr()->isNTerm()) {
      n_term_peak_ptrs.push_back(theo_peak_ptrs[i]);
    }
  }
  return n_term_peak_ptrs;
}

double refinePrecursorAndHeaderShift(ProteoformPtr proteo_ptr,
                                     const ExtendMsPtrVec &ms_three_ptr_vec, 
                                     DiagonalHeaderPtrVec &header_ptrs,
                                     double ppo, double min_mass,
                                     double refine_prec_step_width) {

  double prec_mass = ms_three_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  if (header_ptrs.size() == 0) {
    return prec_mass;
  }
  double tole = ms_three_ptr_vec[0]->getMsHeaderPtr()->getErrorTolerance(ppo);

  int one_side_step_num = 0;
  if (tole > 0) {
    one_side_step_num = (int)std::floor(tole * 2 / refine_prec_step_width);
  }
  int step_num = 1 + 2 * one_side_step_num;

  std::vector<int> counts (step_num);
  int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
  /* for the first size - 1 diagonals, obtain the distribution of errors of
   * c-terminal ions */
  PeakIonPairPtrVec matched_pair_ptrs;
  for(size_t i=0; i< header_ptrs.size() - 1; i++){
    for (size_t j = 0; j < ms_three_ptr_vec.size(); j++) {
      TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
          proteo_ptr, ms_three_ptr_vec[j]->getMsHeaderPtr()->getActivationPtr(),
          header_ptrs, i, min_mass);
      TheoPeakPtrVec c_term_peak_ptrs = getCTermTheoPeakPtrs(theo_peak_ptrs);

      int bgn = header_ptrs[i]->getMatchFirstBpPos()-first_res_pos;
      int end = header_ptrs[i]->getMatchLastBpPos()-first_res_pos;

      PeakIonPairPtrVec pair_ptrs = PeakIonPairFactory::findPairs(ms_three_ptr_vec[j], c_term_peak_ptrs, bgn, end, tole);
      matched_pair_ptrs.insert(matched_pair_ptrs.end(), pair_ptrs.begin(), pair_ptrs.end());
    }
  }

  /* for last diagonal, obtain the distribution of errors of n-term ions */
  for (size_t j = 0; j < ms_three_ptr_vec.size(); j++) {
    TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
        proteo_ptr,
        ms_three_ptr_vec[j]->getMsHeaderPtr()->getActivationPtr(),
        header_ptrs,
        header_ptrs.size() - 1, min_mass);
    TheoPeakPtrVec n_term_peak_ptrs = getNTermTheoPeakPtrs(theo_peak_ptrs);

    int bgn = header_ptrs[header_ptrs.size()-1]->getMatchFirstBpPos()-first_res_pos;
    int end = header_ptrs[header_ptrs.size()-1]->getMatchLastBpPos()-first_res_pos;

    PeakIonPairPtrVec pair_ptrs = PeakIonPairFactory::findPairs(ms_three_ptr_vec[j], n_term_peak_ptrs, bgn, end, tole);
    matched_pair_ptrs.insert(matched_pair_ptrs.end(), pair_ptrs.begin(), pair_ptrs.end());
  }

  for (size_t i = 0; i < matched_pair_ptrs.size(); i++) {
    double diff = matched_pair_ptrs[i]->getTheoPeakPtr()->getPosition() - 
        matched_pair_ptrs[i]->getRealPeakPtr()->getPosition();
    int idx = (int)std::round(diff / refine_prec_step_width) + one_side_step_num; 
    if (idx >= 0 && idx < (int)counts.size()) {
      counts[idx]++;
    }
  }
  /* get median positon */
  int median = matched_pair_ptrs.size()/2;
  int sum = 0;
  int best_pos = one_side_step_num;
  for (size_t i = 0; i < counts.size(); i++) {
    int next_sum = sum + counts[i];
    //LOG_DEBUG("pos " << i << " count " << counts[i]);
    if (sum <= median && median <= next_sum) {
      best_pos = i;
      //LOG_DEBUG("best pos " << best_pos << " count " << counts[best_pos]);
      break;
    }
    sum = next_sum;
  }

  double best_delta = - (best_pos - one_side_step_num) * refine_prec_step_width;
  //LOG_DEBUG("best delta " << best_delta);

  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    header_ptrs[i]->changeOnlyCTermShift(best_delta);
  }
  DiagonalHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size()-1];
  last_header_ptr->changeOnlyNTermShift(best_delta);

  return prec_mass + best_delta;
}

DiagonalHeaderPtrVec refineHeadersBgnEnd(
        ProteoformPtr proteo_ptr,
        const ExtendMsPtrVec &ms_three_ptr_vec,
        const DiagonalHeaderPtrVec& header_ptrs,
        double min_mass){

    DiagonalHeaderPtrVec result_list;
    int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
    int last_res_pos = header_ptrs[header_ptrs.size()-1]->getTruncLastResPos();
    for(size_t i=0; i<header_ptrs.size();i++){
        int bgn = header_ptrs[i]->getMatchFirstBpPos()-first_res_pos;
        int end = header_ptrs[i]->getMatchLastBpPos()-first_res_pos;
        PeakIonPairPtrVec pair_ptrs; 
        for (size_t j = 0; j < ms_three_ptr_vec.size(); j++) {
            TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
                    proteo_ptr,
                    ms_three_ptr_vec[j]->getMsHeaderPtr()->getActivationPtr(),
                    header_ptrs,
                    i, min_mass);
            PeakIonPairPtrVec cur_pair_ptrs = 
                PeakIonPairFactory::findPairs(ms_three_ptr_vec[j], theo_peak_ptrs, bgn, end, 0);
            pair_ptrs.insert(pair_ptrs.end(), cur_pair_ptrs.begin(), cur_pair_ptrs.end());
        }
        if(pair_ptrs.size()<1){
            int pair_size = pair_ptrs.size();
            LOG_TRACE("Empty Segment is found "+StringUtil::convertToString(pair_size));
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

DiagonalHeaderPtrVec2D refineHeadersBgnEnd(
        ProteoformPtr proteo_ptr,
        const ExtendMsPtrVec &ms_three_ptr_vec,
        const DiagonalHeaderPtrVec2D& header_ptrs_2d,
        const DiagonalHeaderPtrVec& header_ptrs_1d,
        double min_mass) {

    DiagonalHeaderPtrVec2D result_list;
    DiagonalHeaderPtr first_header = header_ptrs_1d[0];
    DiagonalHeaderPtr last_header = header_ptrs_1d[header_ptrs_1d.size()-1];
    int first_res_pos = first_header->getTruncFirstResPos();
    int last_res_pos = last_header->getTruncLastResPos();
    int index = 0;
    for(size_t i=0; i<header_ptrs_2d.size();i++){
        DiagonalHeaderPtrVec cur_vec;
        for (size_t j=0; j < header_ptrs_2d[i].size(); j++) {
            int bgn = header_ptrs_2d[i][j]->getMatchFirstBpPos()-first_res_pos;
            int end = header_ptrs_2d[i][j]->getMatchLastBpPos()-first_res_pos;
            PeakIonPairPtrVec pair_ptrs; 
            for (size_t k = 0; k < ms_three_ptr_vec.size(); k++) {
                TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
                        proteo_ptr,
                        ms_three_ptr_vec[k]->getMsHeaderPtr()->getActivationPtr(),
                        header_ptrs_1d,
                        index, min_mass);
                PeakIonPairPtrVec cur_pair_ptrs 
                    = PeakIonPairFactory::findPairs(ms_three_ptr_vec[k], theo_peak_ptrs, bgn, end, 0);
                pair_ptrs.insert(pair_ptrs.end(), cur_pair_ptrs.begin(), cur_pair_ptrs.end());
            }
            index++;
            if(pair_ptrs.size() < 1){
                int pair_size = pair_ptrs.size();
                LOG_TRACE("Empty Segment is found "+ StringUtil::convertToString(pair_size));
                if (header_ptrs_2d[i][j] == first_header) {
                    int new_bgn = first_res_pos;
                    int new_end = first_res_pos;
                    header_ptrs_2d[i][j]->setMatchFirstBpPos(new_bgn);
                    header_ptrs_2d[i][j]->setMatchLastBpPos(new_end);
                    cur_vec.push_back(header_ptrs_2d[i][j]);
                }
                else if (header_ptrs_2d[i][j] == last_header) {
                    int new_bgn = last_res_pos + 1;
                    int new_end = last_res_pos + 1;
                    header_ptrs_2d[i][j]->setMatchFirstBpPos(new_bgn);
                    header_ptrs_2d[i][j]->setMatchLastBpPos(new_end);
                    cur_vec.push_back(header_ptrs_2d[i][j]);
                }
            }
            else{
                int new_bgn = first_res_pos + getNewBgn(pair_ptrs);
                int new_end = first_res_pos + getNewEnd(pair_ptrs);
                header_ptrs_2d[i][j]->setMatchFirstBpPos(new_bgn);
                header_ptrs_2d[i][j]->setMatchLastBpPos(new_end);
                cur_vec.push_back(header_ptrs_2d[i][j]);
            }
        }
        if (cur_vec.size() > 0) {
            result_list.push_back(cur_vec);
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

/*double oldRefinePrecursorAndHeaderShift(ProteoformPtr proteo_ptr,*/
//ExtendMsPtr ms_three_ptr, 
//DiagonalHeaderPtrVec &header_ptrs,
//PtmMngPtr mng_ptr) {
//double prec_mass = ms_three_ptr->getHeaderPtr()->getPrecMonoMass();
//if (header_ptrs.size() == 0) {
//return prec_mass;
//}
//double tole = ms_three_ptr->getHeaderPtr()->getErrorTolerance();
//int one_side_step_num = 0;
//if (tole > 0) {
//one_side_step_num = (int)std::floor(tole/mng_ptr->refine_prec_step_width_);
//}
//int step_num = 1 + 2 * one_side_step_num;
//DiagonalHeaderPtrVec test_ptrs;
//for (size_t i = 0; i < header_ptrs.size(); i++) {
//test_ptrs.push_back(header_ptrs[i]->clone());
//}
//for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
//int middle = (test_ptrs[i]->getMatchLastBpPos() + test_ptrs[i+1]->getMatchFirstBpPos())/2;
//test_ptrs[i]->setMatchLastBpPos(middle);
//test_ptrs[i+1]->setMatchFirstBpPos(middle);
//}
//test_ptrs[0]->setMatchFirstBpPos(test_ptrs[0]->getTruncFirstResPos());
//DiagonalHeaderPtr last_test_ptr = test_ptrs[test_ptrs.size()-1];
//last_test_ptr->setMatchLastBpPos(last_test_ptr->getTruncLastResPos());
////LOG_DEBUG("TEST PTR SIZE" << test_ptrs.size());
//double change = - one_side_step_num * mng_ptr->refine_prec_step_width_;
/* for the first test_ptrs.size() - 1 headers, change C-term shift,
 * for the last header, change N-term shift */
//for (size_t i = 0; i < test_ptrs.size() - 1; i++) {
//test_ptrs[i]->changeOnlyCTermShift(change);
//}
//last_test_ptr->changeOnlyNTermShift(change);
//int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
//double best_score = -1;
//int best_i_start = -1;
//int best_i_end = -1;
//bool continuous = false;
//[>
//for (size_t i = 0; i < ms_three_ptr->size(); i++) {
//std::cout << "ms " << ms_three_ptr->getPeakPtr(i)->getPosition() << std::endl;
//}
//*/
//for (int i = 0; i < step_num; i++) {
////double delta = (i - one_side_step_num) * mng_ptr->refine_prec_step_width_;
//double cur_score = 0;
//for(size_t j=0; j< test_ptrs.size(); j++){
//TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
//proteo_ptr,
//ms_three_ptr->getHeaderPtr()->getActivationPtr(),
//test_ptrs,
//j,
//mng_ptr->prsm_para_ptr_->getSpParaPtr()->getMinMass());

//int bgn = test_ptrs[j]->getMatchFirstBpPos()-first_res_pos;
//int end = test_ptrs[j]->getMatchLastBpPos()-first_res_pos;
//PeakIonPairPtrVec pair_ptrs = findPairs(ms_three_ptr, theo_peak_ptrs, bgn, end, 0 );
//cur_score += pair_ptrs.size();
//}
////LOG_DEBUG("i " << i << " header number "<<  test_ptrs.size() << " delta " << delta << " cur score " << cur_score);
//if (cur_score > best_score) {
//best_score = cur_score;
//best_i_start = i;
//best_i_end = i;
//continuous = true;
//}
//else if (cur_score == best_score) {
//if (continuous) {
//best_i_end = i;
//}
//}
//else {
//continuous = false;
//}

//change = mng_ptr->refine_prec_step_width_;
//for (size_t j = 0; j < test_ptrs.size() - 1; j++) {
//test_ptrs[j]->changeOnlyCTermShift(change);
//}
//last_test_ptr->changeOnlyNTermShift(change);
//}
//for (size_t i = 0; i < header_ptrs.size(); i++) {
//header_ptrs[i]->setMatchFirstBpPos(test_ptrs[i]->getMatchFirstBpPos());
//header_ptrs[i]->setMatchLastBpPos(test_ptrs[i]->getMatchLastBpPos());
//}
//double best_i = (best_i_end + best_i_start) /2; 
//double best_delta = (best_i - one_side_step_num) * mng_ptr->refine_prec_step_width_;
////LOG_DEBUG("best_delta " << best_delta);
//for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
//header_ptrs[i]->changeOnlyCTermShift(best_delta);
//}
//DiagonalHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size()-1];
//last_header_ptr->changeOnlyNTermShift(best_delta);

//return prec_mass + best_delta;
/*}  */

} /* namespace prot */
