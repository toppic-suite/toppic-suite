// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "base/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_top_selector.hpp"

namespace prot {

PrsmTopSelector::PrsmTopSelector(const std::string &db_file_name,
                                 const std::string &spec_file_name,
                                 const std::string &in_file_ext,
                                 const std::string &out_file_ext, 
                                 int n_top): 
    spec_file_name_(spec_file_name), 
    db_file_name_(db_file_name),
    input_file_ext_(in_file_ext),
    output_file_ext_(out_file_ext),
    n_top_(n_top) {
    }

bool containsSameFastaSeq(const PrsmStrPtrVec prsm_ptrs, PrsmStrPtr target_prsm_ptr) {
  for(size_t i=0; i< prsm_ptrs.size();i++){
    if (prsm_ptrs[i]->getSeqName() == target_prsm_ptr->getSeqName()) {
      return true;
    }
  }
  return false;
}

PrsmStrPtrVec getTopPrsms(PrsmStrPtrVec &prsm_str_ptrs, int n_top){
  std::sort(prsm_str_ptrs.begin(),prsm_str_ptrs.end(),PrsmStr::cmpEValueInc);
  int size = prsm_str_ptrs.size();
  int max = size > n_top? n_top:size;
  PrsmStrPtrVec result_ptrs;
  for(int i=0;i<max;i++){
    if(!containsSameFastaSeq(result_ptrs, prsm_str_ptrs[i])){
      result_ptrs.push_back(prsm_str_ptrs[i]);
    }
  }
  return result_ptrs;
}

void PrsmTopSelector::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;
  PrsmReader reader(input_file_name);
  PrsmStrPtr prsm_str_ptr = reader.readOnePrsmStr();

  PrsmXmlWriter writer(base_name +"."+output_file_ext_);
  
  int spec_id = 0;
  while (prsm_str_ptr != nullptr) {
    PrsmStrPtrVec cur_str_ptrs;
    while (prsm_str_ptr != nullptr && prsm_str_ptr->getSpectrumId() == spec_id) {
      cur_str_ptrs.push_back(prsm_str_ptr);
      prsm_str_ptr = reader.readOnePrsmStr();
    }
    PrsmStrPtrVec result_ptrs = getTopPrsms(cur_str_ptrs, n_top_);
    for (size_t i = 0; i < result_ptrs.size(); i++) {
      writer.write(result_ptrs[i]);
    }

    spec_id++;
  }
  
  reader.close();
  writer.close();
}

} /* namespace prot */
