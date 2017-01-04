// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_util.hpp"

namespace prot {

SimplePrsmPtrVec SimplePrsmUtil::getUniqueMatches(SimplePrsmPtrVec match_ptrs) {
  std::sort(match_ptrs.begin(), match_ptrs.end(), SimplePrsm::cmpNameIncScoreDec);
  SimplePrsmPtrVec unique_match_ptrs;
  std::string prev_name = "";
  int prev_id = -1;
  for(size_t i = 0; i < match_ptrs.size(); i++){
    std::string cur_name = match_ptrs[i]->getSeqName();
    int cur_id = match_ptrs[i]->getSpectrumId();
    if (cur_name != prev_name || cur_id != prev_id) {
      unique_match_ptrs.push_back(match_ptrs[i]);
      prev_name = cur_name;
      prev_id = cur_id;
    }
  }

  return unique_match_ptrs;
}


} /* namespace prot */
