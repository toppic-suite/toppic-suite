//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include "common/util/file_util.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_cluster.hpp"

namespace toppic {

namespace prsm_cluster {

std::vector<PrsmStrPtrVec> groupProteins(const PrsmStrPtrVec &prsm_ptrs) {
  // get max shift number
  int max_shift_number = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int cur_shift_number = prsm_ptrs[i]->getUnexpectedPtmNum();
    if (max_shift_number < cur_shift_number) {
      max_shift_number = cur_shift_number;
    }
  }
  // get proteoform groups
  std::vector<PrsmStrPtrVec> proteogroups;
  for (int shift = 0; shift <= max_shift_number; shift++) {
    PrsmStrPtrVec proteo_ptrs;
    for (size_t i = 0; i < prsm_ptrs.size(); i++) {
      if (shift == prsm_ptrs[i]->getUnexpectedPtmNum()) {
        proteo_ptrs.push_back(prsm_ptrs[i]);
      }
    }
    proteogroups.push_back(proteo_ptrs);
  }
  return proteogroups;
}

std::vector<PrsmStrPtrVec> getZeroPtmList(const PrsmStrPtrVec& proteo_ptrs, 
                                          double ppo) {
  std::vector<PrsmStrPtrVec> clusters;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    bool is_found = false;
    for (size_t j = 0; j < clusters.size(); j++) {
      if (PrsmStr::isSameSeqAndMass(proteo_ptrs[i], clusters[j][0], ppo)) {
        clusters[j].push_back(proteo_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmStrPtrVec new_clusters;
      new_clusters.push_back(proteo_ptrs[i]);
      clusters.push_back(new_clusters);
    }
  }
  return clusters;
}

void setProtId(PrsmStrPtrVec& prsm_ptrs) {
  std::vector<PrsmStrPtrVec> proteins;
  std::vector<std::string> protein_names;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    std::string name = prsm_ptrs[i]->getSeqName();
    bool is_found = false;
    for (size_t j = 0; j < protein_names.size(); j++) {
      if (protein_names[j] == name) {
        proteins[j].push_back(prsm_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmStrPtrVec new_protein;
      new_protein.push_back(prsm_ptrs[i]);
      proteins.push_back(new_protein);
      protein_names.push_back(name);
    }
  }

  for (size_t i = 0; i < proteins.size(); i++) {
    for (size_t j = 0; j < proteins[i].size(); j++) {
      proteins[i][j]->setProtId(i);
    }
  }
}

void setClusterId(PrsmStrPtrVec& prsm_ptrs, double ppo) {
  std::vector<PrsmStrPtrVec> proteo_groups = groupProteins(prsm_ptrs);

  // find zero ptm clusters
  std::vector<PrsmStrPtrVec> clusters = getZeroPtmList(proteo_groups[0], ppo);

  for (size_t i = 1; i < proteo_groups.size(); i++) {
    for (size_t j = 0; j < proteo_groups[i].size(); j++) {
      bool is_found = false;
      for (size_t m = 0; m < clusters.size(); m++) {
        if (PrsmStr::isStrictCompatiablePtmSpecies(proteo_groups[i][j], clusters[m][0], ppo)) {
          clusters[m].push_back(proteo_groups[i][j]);
          is_found = true;
          break;
        }
      }
      if (!is_found) {
        PrsmStrPtrVec new_clusters;
        new_clusters.push_back(proteo_groups[i][j]);
        clusters.push_back(new_clusters);
      }
    }
  }

  for (size_t i = 0; i < clusters.size(); i++) {
    for (size_t j = 0; j < clusters[i].size(); j++) {
      clusters[i][j]->setProteoClusterId(i);
    }
  }
}

void process(const std::string &spec_file_name,
             const std::string &input_file_ext,
             const std::string &output_file_ext,
             double ppo)  {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;
  PrsmStrPtrVec prsm_ptrs = prsm_reader_util::readAllPrsmStrs(input_file_name);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpEValueIncProtInc);
  setProtId(prsm_ptrs);
  setClusterId(prsm_ptrs, ppo);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(),
       PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  // output
  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prsm_cluster

}  // namespace toppic


