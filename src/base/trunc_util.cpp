#include "base/logger.hpp"
#include "base/trunc_util.hpp"

namespace prot {

bool TruncUtil::isValidTrunc(TruncPtr trunc_ptr, const ResiduePtrVec & res_ptr_vec) {
  //check if trunc acids match N-terminal acids of the protein 
  int trunc_len = trunc_ptr->getTruncLen();
  if (trunc_len  >= (int)res_ptr_vec.size()) {
    return false;
  }

  ResiduePtrVec trunc_residue_ptr_vec = trunc_ptr->getTruncResiduePtrVec();
  for(int i=0;i<trunc_ptr->getTruncLen();i++){
    if(trunc_residue_ptr_vec[i] != res_ptr_vec[i]){
      return false;
    }
  }
  // check the second letter for NME
  ResiduePtrVec allow_first_remain_residues = trunc_ptr->getAllowFirstRemainResiduePtrs();
  for (size_t i = 0; i < allow_first_remain_residues.size(); i++) {
    if (res_ptr_vec[trunc_len] == allow_first_remain_residues[i]) {
      return true;
    }
  }
  return false;
}

}
