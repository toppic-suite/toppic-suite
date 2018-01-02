//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include "base/logger.hpp"
#include "base/trunc_util.hpp"

namespace prot {

namespace trunc_util {

bool isValidTrunc(TruncPtr trunc_ptr, const ResiduePtrVec & res_ptr_vec) {
  //check if trunc acids match N-terminal acids of the protein 
  int trunc_len = trunc_ptr->getTruncLen();
  if (trunc_len >= (int)res_ptr_vec.size()) {
    return false;
  }

  ResiduePtrVec trunc_residue_ptr_vec = trunc_ptr->getTruncResiduePtrVec();
  for(int i = 0; i < trunc_ptr->getTruncLen(); i++){
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

} // namespace trunc_util

}  // namespace prot
