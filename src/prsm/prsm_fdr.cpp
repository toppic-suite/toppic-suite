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


#include "base/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_fdr.hpp"

namespace prot {
PrsmFdr::PrsmFdr(const std::string &db_file_name,
                 const std::string &spec_file_name,
                 const std::string &input_file_ext,
                 const std::string &output_file_ext): 
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext) {
    }

inline PrsmStrPtr2D getGroups(PrsmStrPtrVec &prsm_ptrs) {
  PrsmStrPtr2D results;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    bool found = false;
    for (size_t j = 0; j < results.size(); j++) {
      if (results[j][0]->getSpeciesId() == prsm_ptrs[i]->getSpeciesId()) {
        found = true;
        results[j].push_back(prsm_ptrs[i]);
        break;
      }
    }
    if (!found) {
      PrsmStrPtrVec new_group;
      new_group.push_back(prsm_ptrs[i]);
      results.push_back(new_group);
    }
  }
  return results;
}

void PrsmFdr::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;

  PrsmStrPtrVec prsm_str_ptrs = PrsmReader::readAllPrsmStrs(input_file_name);

  PrsmStrPtrVec target_ptrs;
  PrsmStrPtrVec decoy_ptrs;
  for(size_t i=0; i< prsm_str_ptrs.size(); i++){
    if (prsm_str_ptrs[i]->getEValue() == 0.0) {
      LOG_ERROR("prot::PRSMFdr zero E value is reported");
    }
    else {
      std::string seq_name  = prsm_str_ptrs[i]->getSeqName();
      //LOG_DEBUG("seq name " << seq_name);
      if(seq_name.find("DECOY_")==0){
        decoy_ptrs.push_back(prsm_str_ptrs[i]);
      }
      else{
        target_ptrs.push_back(prsm_str_ptrs[i]);
      }
    }
  }
  std::sort(target_ptrs.begin(),target_ptrs.end(),PrsmStr::cmpEValueInc);
  std::sort(decoy_ptrs.begin(),decoy_ptrs.end(),PrsmStr::cmpEValueInc);
  
  computeFdr(target_ptrs,decoy_ptrs);

  PrsmStrPtr2D target_proteoforms = getGroups(target_ptrs);
  PrsmStrPtr2D decoy_proteoforms = getGroups(decoy_ptrs);

  computeProteoformFdr(target_proteoforms, decoy_proteoforms);

  std::string output_file_name = base_name+"."+output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  std::sort(target_ptrs.begin(),target_ptrs.end(),PrsmStr::cmpSpectrumIdInc);
  writer.writeVector(target_ptrs);
  writer.close();
}

void PrsmFdr::computeFdr(PrsmStrPtrVec &target_ptrs,PrsmStrPtrVec &decoy_ptrs){
  for(size_t i=0; i<target_ptrs.size(); i++){
    int n_target=i+1;
    double target_evalue = target_ptrs[i]->getEValue();
    int n_decoy = 0;
    for(size_t j=0; j<decoy_ptrs.size(); j++){
      if(decoy_ptrs[j]->getEValue() <= target_evalue){
        n_decoy++;
      }
      else{
        break;
      }
    }
    double fdr = (double)n_decoy/(double)n_target;
    if(fdr>1){
      fdr=1.0;
    }
    target_ptrs[i]->setFdr(fdr);
  }
}

void PrsmFdr::computeProteoformFdr(PrsmStrPtr2D &target_proteoforms,
                                   PrsmStrPtr2D &decoy_proteoforms) {
  for(size_t i=0; i<target_proteoforms.size(); i++){
    int n_target=i+1;
    double target_evalue = target_proteoforms[i][0]->getEValue();
    int n_decoy = 0;
    for(size_t j=0; j<decoy_proteoforms.size(); j++){
      if(decoy_proteoforms[j][0]->getEValue() <= target_evalue){
        n_decoy++;
      }
      else{
        break;
      }
    }
    double fdr = (double)n_decoy/(double)n_target;
    if(fdr>1){
      fdr=1.0;
    }
    for (size_t k = 0; k < target_proteoforms[i].size(); k++) {
      target_proteoforms[i][k]->setProteoformFdr(fdr);
    }
  }
}

} /* namespace prot */
