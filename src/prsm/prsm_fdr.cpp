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


#include "util/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_fdr.hpp"

namespace toppic {

inline PrsmStrPtrVec2D getGroups(PrsmStrPtrVec &prsm_ptrs) {
  PrsmStrPtrVec2D results;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    bool found = false;
    for (size_t j = 0; j < results.size(); j++) {
      if (results[j][0]->getClusterId() == prsm_ptrs[i]->getClusterId()) {
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
  std::string base_name = file_util::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;

  PrsmStrPtrVec prsm_str_ptrs = PrsmReader::readAllPrsmStrs(input_file_name);

  PrsmStrPtrVec target_ptrs;
  PrsmStrPtrVec decoy_ptrs;
  for(size_t i = 0; i < prsm_str_ptrs.size(); i++){
    if (prsm_str_ptrs[i]->getEValue() == 0.0) {
      LOG_ERROR("toppic::PRSMFdr zero E value is reported");
    } else {
      std::string seq_name  = prsm_str_ptrs[i]->getSeqName();
      //LOG_DEBUG("seq name " << seq_name);
      if(seq_name.find("DECOY_")==0){
        decoy_ptrs.push_back(prsm_str_ptrs[i]);
      } else{
        target_ptrs.push_back(prsm_str_ptrs[i]);
      }
    }
  }
  std::sort(target_ptrs.begin(),target_ptrs.end(),PrsmStr::cmpEValueInc);
  std::sort(decoy_ptrs.begin(),decoy_ptrs.end(),PrsmStr::cmpEValueInc);

  computeFdr(target_ptrs,decoy_ptrs);

  PrsmStrPtrVec2D target_proteoforms = getGroups(target_ptrs);
  PrsmStrPtrVec2D decoy_proteoforms = getGroups(decoy_ptrs);

  computeProteoformFdr(target_proteoforms, decoy_proteoforms);

  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  std::sort(target_ptrs.begin(), target_ptrs.end(), PrsmStr::cmpSpectrumIdInc);
  writer.writeVector(target_ptrs);
  writer.close();
}

void PrsmFdr::computeFdr(PrsmStrPtrVec &target_ptrs, PrsmStrPtrVec &decoy_ptrs){
  for(size_t i = 0; i < target_ptrs.size(); i++){
    int n_target = i + 1;
    double target_evalue = target_ptrs[i]->getEValue();
    int n_decoy = 0;
    for(size_t j = 0; j < decoy_ptrs.size(); j++){
      if(decoy_ptrs[j]->getEValue() <= target_evalue){
        n_decoy++;
      } else{
        break;
      }
    }
    double fdr = static_cast<double>(n_decoy) / static_cast<double>(n_target);
    if (fdr > 1) {
      fdr = 1.0;
    }
    target_ptrs[i]->setFdr(fdr);
  }
}

void PrsmFdr::computeProteoformFdr(PrsmStrPtrVec2D &target_proteoforms,
                                   PrsmStrPtrVec2D &decoy_proteoforms) {
  for(size_t i = 0; i < target_proteoforms.size(); i++) {
    int n_target= i + 1;
    double target_evalue = target_proteoforms[i][0]->getEValue();
    int n_decoy = 0;
    for (size_t j = 0; j < decoy_proteoforms.size(); j++) {
      if (decoy_proteoforms[j][0]->getEValue() <= target_evalue) {
        n_decoy++;
      } else {
        break;
      }
    }
    double fdr = static_cast<double>(n_decoy) / static_cast<double>(n_target);
    if(fdr > 1){
      fdr = 1.0;
    }
    for (size_t k = 0; k < target_proteoforms[i].size(); k++) {
      target_proteoforms[i][k]->setProteoformFdr(fdr);
    }
  }
}

}  // namespace toppic
