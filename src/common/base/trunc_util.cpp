//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include "common/util/logger.hpp"
#include "common/base/trunc_util.hpp"

namespace toppic {

namespace trunc_util {

bool isValidTrunc(TruncPtr trunc_ptr, const ResiduePtrVec & res_ptr_vec) {
  // check if trunc acids match N-terminal acids of the protein 
  int trunc_len = trunc_ptr->getTruncLen();
  if (trunc_len >= static_cast<int>(res_ptr_vec.size())) {
    return false;
  }

  ResiduePtrVec trunc_residue_ptr_vec = trunc_ptr->getTruncResiduePtrVec();
  for (int i = 0; i < trunc_ptr->getTruncLen(); i++) {
    // check amino acid match only
    if (trunc_residue_ptr_vec[i]->getAminoAcidPtr() != res_ptr_vec[i]->getAminoAcidPtr()) {
      return false;
    }
  }
  // check the second letter for NME
  ResiduePtrVec allow_first_remain_residues = trunc_ptr->getAllowFirstRemainResiduePtrs();
  for (size_t i = 0; i < allow_first_remain_residues.size(); i++) {
    // check amino acid match only
    if (res_ptr_vec[trunc_len]->getAminoAcidPtr() 
        == allow_first_remain_residues[i]->getAminoAcidPtr()) {
      return true;
    }
  }
  return false;
}

}  // namespace trunc_util

}  // namespace toppic
