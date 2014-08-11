#include <algorithm>
#include <iostream>
#include "filterdiagonal/ptm_fast_filter_hi_mem.hpp"

namespace prot {

PtmFastFilterHiMem::PtmFastFilterHiMem(const ProteoformPtrVec &proteo_ptrs,
                                       PtmFastFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  index_ptr_ = CompShiftHiMemPtr(new CompShiftHiMem(proteo_ptrs, mng_ptr));
}

SimplePrsmPtrVec PtmFastFilterHiMem::getBestMatch(PrmMsPtr ms_ptr){
  SimplePrsmPtrVec2D match_ptrs = compute(ms_ptr);
  SimplePrsmPtrVec unique_match_ptrs = sort(match_ptrs);
  size_t num = mng_ptr_->ptm_fast_filter_result_num_;
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

SimplePrsmPtrVec2D PtmFastFilterHiMem::compute(PrmMsPtr ms_ptr){
  std::pair<std::vector<int>, std::vector<int>> mass_errors 
      = getIntMassErrorList(ms_ptr, mng_ptr_->ptm_fast_filter_scale_,true,false);
  SimplePrsmPtrVec2D match_ptrs;
  for(size_t i=0;i<mass_errors.first.size();i++){
    std::vector<std::pair<int,int>> results 
        =index_ptr_->compConvolution(mass_errors.first, mass_errors.second, i,
                                 mng_ptr_->ptm_fast_filter_result_num_);
    SimplePrsmPtrVec temp_match_ptrs;
    for(size_t j =0;j <results.size();j++){
      temp_match_ptrs.push_back(
          SimplePrsmPtr(new SimplePrsm(ms_ptr->getHeaderPtr(),
                                       proteo_ptrs_[results[j].first],
                                       results[j].second)));
    }
    match_ptrs.push_back(temp_match_ptrs);
  }
  return match_ptrs;
}

SimplePrsmPtrVec PtmFastFilterHiMem::sort(const SimplePrsmPtrVec2D &match_ptrs) {
  SimplePrsmPtrVec sorted_match_ptrs;

  for(size_t i=0;i<match_ptrs.size();i++){
    for(size_t j =0;j< match_ptrs[i].size();j++){
      sorted_match_ptrs.push_back(match_ptrs[i][j]);
    }
  }

  std::sort(sorted_match_ptrs.begin(),sorted_match_ptrs.end(),simplePrsmDown);

  SimplePrsmPtrVec unique_match_ptrs;
  for(size_t i=0;i< sorted_match_ptrs.size();i++){
    bool found = false;
    std::string seq_name = sorted_match_ptrs[i]->getSeqName();
    for(size_t j=0;j<unique_match_ptrs.size();j++){
      if(seq_name == unique_match_ptrs[j]->getSeqName()){
        found=true;
        break;
      }
    }
    if(!found){
      unique_match_ptrs.push_back(sorted_match_ptrs[i]);
    }
  }
  return unique_match_ptrs;
}

} /* namespace prot */
