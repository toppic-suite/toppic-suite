//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_simple_cluster.hpp"

namespace toppic {

PrsmSimpleCluster::PrsmSimpleCluster(const std::string &db_file_name,
                                     const std::string &spec_file_name,
                                     const std::string &input_file_ext,
                                     const std::string &output_file_ext,
                                     double tolerance):
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext),
    tolerance_(tolerance) {}

void PrsmSimpleCluster::setProtId(PrsmStrPtrVec& prsm_ptrs) {
  std::vector<PrsmStrPtrVec> proteins;
  std::vector<std::string> protein_names;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    std::cout << "i" << i << " out of " << prsm_ptrs.size() << std::endl;
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

void PrsmSimpleCluster::setClusterId(PrsmStrPtrVec& prsm_ptrs, double tolerance) {
  std::vector<PrsmStrPtrVec> clusters;
  for (size_t j = 0; j < prsm_ptrs.size(); j++) {
    std::cout << "j" << j << " out of " << prsm_ptrs.size() << std::endl;
    bool is_found = false;
    for (size_t m = 0; m < clusters.size(); m++) {
      if (PrsmStr::isSimpleMatch(prsm_ptrs[j], clusters[m][0], tolerance)) {
        clusters[m].push_back(prsm_ptrs[j]);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmStrPtrVec new_clusters;
      new_clusters.push_back(prsm_ptrs[j]);
      clusters.push_back(new_clusters);
    }
  }

  for (size_t i = 0; i < clusters.size(); i++) {
    for (size_t j = 0; j < clusters[i].size(); j++) {
      clusters[i][j]->setClusterId(i);
    }
  }
}

void PrsmSimpleCluster::process() {
  std::string base_name = file_util::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;
  LOG_DEBUG("Reading prsm strings started");
  PrsmStrPtrVec prsm_ptrs = PrsmReader::readAllPrsmStrs(input_file_name);
  LOG_DEBUG("Reading prsm strings finished");
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpEValueInc);
  setProtId(prsm_ptrs);
  setClusterId(prsm_ptrs, tolerance_);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpSpectrumIdIncPrecursorIdInc);
  // output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace toppic


