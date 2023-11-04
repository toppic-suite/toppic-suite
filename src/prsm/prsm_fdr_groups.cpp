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
#include "prsm/prsm_fdr_groups.hpp"

namespace toppic {

namespace prsm_fdr_groups {

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
  PrsmStrPtrVec all_ptrs1;
  PrsmStrPtrVec all_ptrs2;
  PrsmStrPtrVec all_ptrs3;
  PrsmStrPtrVec all_ptrs4;

  //Group 1: No mass shift and 0 variable PTMS
  PrsmStrPtrVec target_ptrs1;
  PrsmStrPtrVec decoy_ptrs1;
  //Group2: No mass shift and >0 variable PTMs
  PrsmStrPtrVec target_ptrs2;
  PrsmStrPtrVec decoy_ptrs2;
  //Group3: 1 mass shift
  PrsmStrPtrVec target_ptrs3;
  PrsmStrPtrVec decoy_ptrs3;
  //Group4: 2 mass shift
  PrsmStrPtrVec target_ptrs4;
  PrsmStrPtrVec decoy_ptrs4;
  for(size_t i = 0; i < prsm_str_ptrs.size(); i++){
    if (prsm_str_ptrs[i]->getEValue() == 0.0) {
      LOG_ERROR("toppic::PRSMFdr zero E value is reported")
    } else {
      std::string seq_name  = prsm_str_ptrs[i]->getSeqName();
      //LOG_DEBUG("seq name " << seq_name);
      if(seq_name.find("DECOY_")==0){
          //0 mass shift
          if(prsm_str_ptrs[i]->getUnexpectedPtmNum() == 0) {
              //0 variable PTMs
              if (prsm_str_ptrs[i]->getVariablePtmNum() == 0) {
                  decoy_ptrs1.push_back(prsm_str_ptrs[i]);
              }
              //>0 variable PTMs
              else{
                  decoy_ptrs2.push_back(prsm_str_ptrs[i]);
              }
          }
          //1 mass shift
          else if (prsm_str_ptrs[i]->getUnexpectedPtmNum() == 1) {
              decoy_ptrs3.push_back(prsm_str_ptrs[i]);
          }
          //2 mass shift
          else{
              decoy_ptrs4.push_back(prsm_str_ptrs[i]);
          }
      } else{
          //0 mass shift
          if(prsm_str_ptrs[i]->getUnexpectedPtmNum() == 0) {
              //0 variable PTMs
              if (prsm_str_ptrs[i]->getVariablePtmNum() == 0) {
                  target_ptrs1.push_back(prsm_str_ptrs[i]);
              }
              //>0 variable PTMs
              else{
                  target_ptrs2.push_back(prsm_str_ptrs[i]);
              }
          }
          //1 mass shift
          else if (prsm_str_ptrs[i]->getUnexpectedPtmNum() == 1) {
              target_ptrs3.push_back(prsm_str_ptrs[i]);
          }
          //2 mass shift
          else{
              target_ptrs4.push_back(prsm_str_ptrs[i]);
          }
      }
    }
  }
  std::sort(target_ptrs1.begin(),target_ptrs1.end(),PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_ptrs1.begin(),decoy_ptrs1.end(),PrsmStr::cmpEValueIncProtInc);
  std::sort(target_ptrs2.begin(),target_ptrs2.end(),PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_ptrs2.begin(),decoy_ptrs2.end(),PrsmStr::cmpEValueIncProtInc);
  std::sort(target_ptrs3.begin(),target_ptrs3.end(),PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_ptrs3.begin(),decoy_ptrs3.end(),PrsmStr::cmpEValueIncProtInc);
  std::sort(target_ptrs4.begin(),target_ptrs4.end(),PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_ptrs4.begin(),decoy_ptrs4.end(),PrsmStr::cmpEValueIncProtInc);

  computeFdr(target_ptrs1,decoy_ptrs1);
  computeFdr(target_ptrs2,decoy_ptrs2);
  computeFdr(target_ptrs3,decoy_ptrs3);
  computeFdr(target_ptrs4,decoy_ptrs4);


  PrsmStrPtrVec2D target1_proteoforms = getGroups(target_ptrs1);
  PrsmStrPtrVec2D decoy1_proteoforms = getGroups(decoy_ptrs1);
  PrsmStrPtrVec2D target2_proteoforms = getGroups(target_ptrs2);
  PrsmStrPtrVec2D decoy2_proteoforms = getGroups(decoy_ptrs2);
  PrsmStrPtrVec2D target3_proteoforms = getGroups(target_ptrs3);
  PrsmStrPtrVec2D decoy3_proteoforms = getGroups(decoy_ptrs3);
  PrsmStrPtrVec2D target4_proteoforms = getGroups(target_ptrs4);
  PrsmStrPtrVec2D decoy4_proteoforms = getGroups(decoy_ptrs4);

  computeProteoformFdr(target1_proteoforms, decoy1_proteoforms);
  computeProteoformFdr(target2_proteoforms, decoy2_proteoforms);
  computeProteoformFdr(target3_proteoforms, decoy3_proteoforms);
  computeProteoformFdr(target4_proteoforms, decoy4_proteoforms);


  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);

  if (keep_decoy_results == "true") {
    //concat target_ptrs and decoy_ptrs;
    all_ptrs1.insert(all_ptrs1.begin(), target_ptrs1.begin(), target_ptrs1.end());
    all_ptrs1.insert(all_ptrs1.end(), decoy_ptrs1.begin(), decoy_ptrs1.end());
    std::sort(all_ptrs1.begin(), all_ptrs1.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    std::cout << "Group 1 target_ptrs: " << target_ptrs1.size() << ", Group 1 decoy_ptrs: " << decoy_ptrs1.size() << std::endl;
    writer.writeVector(all_ptrs1);

    all_ptrs2.insert(all_ptrs2.begin(), target_ptrs2.begin(), target_ptrs2.end());
    all_ptrs2.insert(all_ptrs2.end(), decoy_ptrs2.begin(), decoy_ptrs2.end());
    std::sort(all_ptrs2.begin(), all_ptrs2.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    std::cout << "Group 2 target_ptrs: " << target_ptrs2.size() << ", Group 2 decoy_ptrs: " << decoy_ptrs2.size() << std::endl;
    writer.writeVector(all_ptrs2);

    all_ptrs3.insert(all_ptrs3.begin(), target_ptrs3.begin(), target_ptrs3.end());
    all_ptrs3.insert(all_ptrs3.end(), decoy_ptrs3.begin(), decoy_ptrs3.end());
    std::sort(all_ptrs3.begin(), all_ptrs3.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    std::cout << "Goup 3 target_ptrs: " << target_ptrs3.size() << ", Group 3 decoy_ptrs: " << decoy_ptrs3.size() << std::endl;
    writer.writeVector(all_ptrs3);

    all_ptrs4.insert(all_ptrs4.begin(), target_ptrs4.begin(), target_ptrs4.end());
    all_ptrs4.insert(all_ptrs4.end(), decoy_ptrs4.begin(), decoy_ptrs4.end());
    std::sort(all_ptrs4.begin(), all_ptrs4.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    std::cout << "Goup 4 target_ptrs: " << target_ptrs4.size() << ", Group 4 decoy_ptrs: " << decoy_ptrs4.size() << std::endl;
    writer.writeVector(all_ptrs4);
  }
  else {
    std::sort(target_ptrs1.begin(), target_ptrs1.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    writer.writeVector(target_ptrs1);

    std::sort(target_ptrs2.begin(), target_ptrs2.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    writer.writeVector(target_ptrs2);

    std::sort(target_ptrs3.begin(), target_ptrs3.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    writer.writeVector(target_ptrs3);

    std::sort(target_ptrs4.begin(), target_ptrs4.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    writer.writeVector(target_ptrs4);
  }
  writer.close();
}

}

}  // namespace toppic
