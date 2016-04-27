#include <algorithm>
#include <iostream>
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_match_factory.hpp"
#include "zeroptmfilter/mass_match_util.hpp"
#include "diagfilter/mass_diag_filter.hpp"

namespace prot {

MassDiagFilter::MassDiagFilter(const ProteoformPtrVec &proteo_ptrs,
                               DiagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  bool rev = false;
  index_ptr_ = MassMatchFactory::getMassMatchPtr(proteo_ptrs, mng_ptr->max_proteoform_mass_, 
                                                 mng_ptr->filter_scale_, rev);
}

SimplePrsmPtrVec MassDiagFilter::getBestMatch(const PrmMsPtrVec &ms_ptr_vec){
  SimplePrsmPtrVec match_ptrs = compute(ms_ptr_vec);
  SimplePrsmPtrVec unique_match_ptrs = SimplePrsmUtil::getUniqueMatches(match_ptrs);
  std::sort(unique_match_ptrs.begin(), unique_match_ptrs.end(), SimplePrsm::cmpScoreDec);
  size_t num = mng_ptr_->filter_result_num_;
  if(num > unique_match_ptrs.size()){
    num = unique_match_ptrs.size();
  }
  SimplePrsmPtrVec result_ptrs;
  for(size_t i=0;i<num;i++){
    SimplePrsmPtr match_ptr = unique_match_ptrs[i];
    if(match_ptr->getScore() > 0.0){
      result_ptrs.push_back(match_ptr);
    }
    else{
      break;
    }
  }
  return result_ptrs;
}

inline SimplePrsmPtrVec MassDiagFilter::compute(const PrmMsPtrVec &ms_ptr_vec){
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> mass_errors 
      = PrmMs::getIntMassErrorList(ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_,true,false);
  //LOG_DEBUG("mass error size " << mass_errors.size() << " filter result number " << mng_ptr_->filter_result_num_);
  SimplePrsmPtrVec match_ptrs;
  int row_num = index_ptr_->getRowNum();
  std::vector<int> proteo_row_begins = index_ptr_->getProteoRowBegins();
  std::vector<int> proteo_row_ends = index_ptr_->getProteoRowEnds();
  int threshold = 4;
  for(size_t i=0;i<mass_errors.size();i++){
    std::vector<short> scores(row_num, 0);
    index_ptr_->compScores(mass_errors, i, -mass_errors[i].first, scores);
    FilterProteinPtrVec results 
        = MassMatchUtil::findTopProteins(scores, proteo_row_begins, proteo_row_ends, threshold,
                                         mng_ptr_->filter_result_num_);
    //LOG_DEBUG("result size " << results.size());
    for(size_t j =0;j <results.size();j++){
      int id = results[j]->getProteinId();
      int score = results[j]->getScore();
      match_ptrs.push_back(
          SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(),
                                       ms_ptr_vec.size(),
                                       proteo_ptrs_[id], score)));
    }
  }
  return match_ptrs;
}

} /* namespace prot */
