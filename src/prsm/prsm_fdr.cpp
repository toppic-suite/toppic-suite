//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_fdr.hpp"

namespace toppic {

namespace prsm_fdr {

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

double computeFdr(int n_decoy, int n_target) {
  double fdr = static_cast<double>(n_decoy) / static_cast<double>(n_target);
  if(fdr > 1){
    fdr = 1.0;
  }
  return fdr;
}


void computeFdr(PrsmStrPtrVec &target_ptrs, PrsmStrPtrVec &decoy_ptrs){
  int n_decoy = 0;
  for(size_t i = 0; i < target_ptrs.size(); i++){
    int n_target = i + 1;
    double target_evalue = target_ptrs[i]->getEValue();
    for(size_t j = n_decoy; j < decoy_ptrs.size(); j++){
      if(decoy_ptrs[j]->getEValue() <= target_evalue){
        n_decoy++;
        double fdr = computeFdr(n_decoy, n_target);
        decoy_ptrs[j]->setFdr(fdr);
      } else{
        break;
      }
    }
    double fdr = computeFdr(n_decoy, n_target);
    target_ptrs[i]->setFdr(fdr);
  }
}

void computeProteoformFdr(PrsmStrPtrVec2D &target_proteoforms,
                          PrsmStrPtrVec2D &decoy_proteoforms) {
  int n_decoy = 0;
  for(size_t i = 0; i < target_proteoforms.size(); i++) {
    int n_target= i + 1;
    double target_evalue = target_proteoforms[i][0]->getEValue();
    for (size_t j = n_decoy; j < decoy_proteoforms.size(); j++) {
      if (decoy_proteoforms[j][0]->getEValue() <= target_evalue) {
        n_decoy++;
        double fdr = computeFdr(n_decoy, n_target);
        for (size_t k = 0; k < decoy_proteoforms[j].size(); k++) {
          decoy_proteoforms[j][k]->setProteoformFdr(fdr);
        }
      } else {
        break;
      }
    }
    double fdr = computeFdr(n_decoy, n_target);
    for (size_t k = 0; k < target_proteoforms[i].size(); k++) {
      target_proteoforms[i][k]->setProteoformFdr(fdr);
    }
  }
}


void process(const std::string &spec_file_name,
             const std::string &input_file_ext,
             const std::string &output_file_ext, 
             std::string keep_decoy_results) { 
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;

  PrsmStrPtrVec prsm_str_ptrs = prsm_reader_util::readAllPrsmStrs(input_file_name);
  PrsmStrPtrVec all_ptrs;
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

  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);

  if (keep_decoy_results == "true") {
    //concat target_ptrs and decoy_ptrs;
    all_ptrs.insert(all_ptrs.begin(), target_ptrs.begin(), target_ptrs.end());
    all_ptrs.insert(all_ptrs.end(), decoy_ptrs.begin(), decoy_ptrs.end());
    std::sort(all_ptrs.begin(), all_ptrs.end(), PrsmStr::cmpSpectrumIdInc);
    std::cout << "target_ptrs: " << target_ptrs.size() << ", decoy_ptrs: " << decoy_ptrs.size() << std::endl;
    writer.writeVector(all_ptrs);
  }
  else {
    std::sort(target_ptrs.begin(), target_ptrs.end(), PrsmStr::cmpSpectrumIdInc);
    writer.writeVector(target_ptrs);
  }
  writer.close();
}

}

}  // namespace toppic
