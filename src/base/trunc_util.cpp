#include "base/logger.hpp"
#include "base/trunc_util.hpp"

namespace prot {

bool TruncUtil::isSameTrunc(TruncPtr trunc_ptr, const ResiduePtrVec& res_ptr_vec, int len) {
  if(trunc_ptr->getTruncLen() != len){
    return false;
  }
  AcidPtrVec acid_ptr_vec = trunc_ptr->getAcidPtrVec();
  for(int i=0;i<trunc_ptr->getTruncLen();i++){
    if(acid_ptr_vec[i] != res_ptr_vec[i]->getAcidPtr()){
      return false;
    }
  }
  return true;
}

bool TruncUtil::isValidTrunc(TruncPtr trunc_ptr, const ResiduePtrVec & res_ptr_vec) {
  //check if trunc acids match N-terminal acids of the protein 
  bool result = true;
  int trunc_len = trunc_ptr->getTruncLen();
  if (trunc_len >= (int)res_ptr_vec.size()) {
    result = false; 
  }
  else {
    result = isSameTrunc(trunc_ptr, res_ptr_vec, trunc_len);
  }
  //LOG_DEBUG("Valid trunc " << result << " trunc len " << trunc_len_ 
  //          << " seq len " << res_seq_ptr->getLen());
  return result;
}

}
