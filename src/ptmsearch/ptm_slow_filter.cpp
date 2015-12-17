#include "base/prot_mod_base.hpp"
#include "base/proteoform_factory.hpp"
#include "ptmsearch/ptm_slow_filter.hpp"

namespace prot {

PtmSlowFilter::PtmSlowFilter(
    SpectrumSetPtr spectrum_set_ptr,
    SimplePrsmPtrVec simple_prsm_ptrs,
    CompShiftLowMemPtr comp_shift_ptr,
    PtmMngPtr mng_ptr){

  std::string db_file_name = mng_ptr->prsm_para_ptr_->getSearchDbFileName();
  FastaIndexReaderPtr reader_ptr(new FastaIndexReader(db_file_name));
  ModPtrVec fix_mod_list = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
  ProtModPtrVec prot_mod_ptr_vec = mng_ptr->prsm_para_ptr_->getProtModPtrVec();
  ProteoformPtrVec proteoform_ptr_vec;
  for (size_t i = 0; i < simple_prsm_ptrs.size(); i++) {
    std::string seq_name = simple_prsm_ptrs[i]->getSeqName();
    std::string seq_desc = simple_prsm_ptrs[i]->getSeqDesc();
    ProteoformPtr proteo_ptr = ProteoformFactory::readFastaToProteoformPtr(
        reader_ptr, seq_name, seq_desc, fix_mod_list);
    ProteoformPtrVec mod_form_ptr_vec = ProteoformFactory::geneProtModProteoform(
        proteo_ptr, prot_mod_ptr_vec);
    for (size_t j = 0; j < mod_form_ptr_vec.size(); j++) {
      PtmSlowMatchPtr ptm_slow_match_ptr(new PtmSlowMatch(mod_form_ptr_vec[j],spectrum_set_ptr,
                                                          comp_shift_ptr,mng_ptr));
      complete_prefix_slow_match_ptrs_.push_back(ptm_slow_match_ptr);
    }
  }

  // init suffix_internal_slow_prsm_ptrs
  for(size_t i=0; i<complete_prefix_slow_match_ptrs_.size();i++){
    ProtModPtr prot_mod_ptr = complete_prefix_slow_match_ptrs_[i]->getProteoform()->getProtModPtr();
    if (prot_mod_ptr == ProtModBase::getProtModPtr_NONE()) {
      suffix_internal_slow_match_ptrs_.push_back(complete_prefix_slow_match_ptrs_[i]);
    }
  }

  // compute complete and prefix prsms 
  for(size_t i=0; i<complete_prefix_slow_match_ptrs_.size();i++){
    PrsmPtrVec comp_ptrs;
    //LOG_DEBUG("Compute complete prsm " << i);
    complete_prefix_slow_match_ptrs_[i]->compute(AlignType::COMPLETE, comp_ptrs);
    complete_prsm_2d_ptrs_.push_back(comp_ptrs);
    PrsmPtrVec prefix_ptrs;
    //LOG_DEBUG("Compute prefix prsm " << i);
    complete_prefix_slow_match_ptrs_[i]->compute(AlignType::PREFIX, prefix_ptrs);
    prefix_prsm_2d_ptrs_.push_back(prefix_ptrs);
    //LOG_DEBUG("compute prefi completed");
  }
  //LOG_DEBUG("complete prefix completed");

  // compute suffix and internal prsms 
  for(size_t i=0; i< suffix_internal_slow_match_ptrs_.size();i++){
    PrsmPtrVec suffix_ptrs;
    suffix_internal_slow_match_ptrs_[i]->compute(AlignType::SUFFIX, suffix_ptrs);
    suffix_prsm_2d_ptrs_.push_back(suffix_ptrs);
    PrsmPtrVec internal_ptrs;
    suffix_internal_slow_match_ptrs_[i]->compute(AlignType::INTERNAL, internal_ptrs);
    internal_prsm_2d_ptrs_.push_back(internal_ptrs);
  }
  //LOG_DEBUG("suffix internal completed");
}

PrsmPtrVec PtmSlowFilter::getPrsms(int shift_num, AlignTypePtr type_ptr){
  PrsmPtrVec prsm_ptrs;
  if (type_ptr == AlignType::COMPLETE) {
    for (size_t i = 0; i < complete_prsm_2d_ptrs_.size(); i++) {
      if (complete_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(complete_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  else if (type_ptr == AlignType::PREFIX) {
    for (size_t i = 0; i < prefix_prsm_2d_ptrs_.size(); i++) {
      if (prefix_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(prefix_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  else if (type_ptr == AlignType::SUFFIX) {
    for (size_t i = 0; i < suffix_prsm_2d_ptrs_.size(); i++) {
      if (suffix_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(suffix_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  else if (type_ptr == AlignType::INTERNAL) {
    for (size_t i = 0; i < internal_prsm_2d_ptrs_.size(); i++) {
      if (internal_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(internal_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  return prsm_ptrs;
}

} /* namespace prot */
