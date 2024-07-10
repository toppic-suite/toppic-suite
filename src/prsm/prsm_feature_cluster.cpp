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
#include <iomanip>
#include <iostream>
#include <sstream> 

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_feature_cluster.hpp"

namespace toppic {

namespace prsm_feature_cluster {

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

void setProteoClusterId(PrsmStrPtrVec& prsm_ptrs,
                        double prec_error_tole) {
  std::vector<PrsmStrPtrVec> clusters;
  int prsm_count = prsm_ptrs.size();
  // first round, cluster proteoforms using feature id
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    bool is_found = false;
    PrsmStrPtr cur_ptr = prsm_ptrs[i];
    for (size_t j = 0; j < clusters.size(); j++) {
      PrsmStrPtr ref_ptr = clusters[j][0];
      // if the same feature id
      if (cur_ptr->getFracFeatureId() == ref_ptr->getFracFeatureId()) {
        clusters[j].push_back(cur_ptr);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmStrPtrVec new_clusters;
      new_clusters.push_back(cur_ptr);
      clusters.push_back(new_clusters);
    }
    double perc = (i + 1) * 100.0 / prsm_count;
    std::cout << std::flush << "Finding PrSM clusters - processing " 
        << std::setprecision(3) << perc << "%. \r";
  }
  // second round, merge feature clusters
  std::vector<PrsmStrPtrVec> merged_clusters;
  for (size_t i = 0; i < clusters.size(); i++) {
    bool is_found = false;
    PrsmStrPtr cur_ptr = clusters[i][0];
    for (size_t j = 0; j < merged_clusters.size(); j++) {
      PrsmStrPtr ref_ptr = merged_clusters[j][0];
      // if the same protein and similar mass
      if (cur_ptr->getProtId() == ref_ptr->getProtId()) {
        if (std::abs(cur_ptr->getOriPrecMass() - ref_ptr->getOriPrecMass()) 
            <= prec_error_tole) {
          merged_clusters[j].insert(merged_clusters[j].end(),
                                    clusters[i].begin(), 
                                    clusters[i].end());
          //LOG_DEBUG("Proteoform merging by mass!");
          is_found = true;
          break;
        }
      } 
      else {
        // if protein identifications are different, but the protein sequences
        // are the same and the proteoform masses are similar, the two
        // proteoforms are treated as one. 
        if (cur_ptr->getProteoformDbSeq() == ref_ptr->getProteoformDbSeq()
            && std::abs(cur_ptr->getOriPrecMass() - ref_ptr->getOriPrecMass()) 
            <= prec_error_tole) {
          merged_clusters[j].insert(merged_clusters[j].end(),
                                    clusters[i].begin(), 
                                    clusters[i].end());
          //LOG_DEBUG("Proteoform merging by sequence!");
          is_found = true;
          break;
        }
      }
    }
    if (!is_found) {
      merged_clusters.push_back(clusters[i]);
    }
  }
  std::cout << std::endl;
  for (size_t i = 0; i < merged_clusters.size(); i++) {
    double inte = prsm_util::compClusterInte(merged_clusters[i]);
    for (size_t j = 0; j < merged_clusters[i].size(); j++) {
      merged_clusters[i][j]->setProteoClusterId(i);
      merged_clusters[i][j]->setProteoInte(inte);
    }
  }
}

void process(const std::string &spec_file_name,
             const std::string &input_file_ext,
             const std::string &output_file_ext,
             double prec_error_tole) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;
  PrsmStrPtrVec prsm_ptrs = prsm_reader_util::readAllPrsmStrsMatchSeq(input_file_name);

  std::string feature_file_name = base_name + ".feature";
  prsm_util::addFeatureInfoToPrsms(prsm_ptrs, feature_file_name);
  // remove prsms without feature
  PrsmStrPtrVec filtered_prsm_ptrs;
  prsm_util::removePrsmsWithoutFeature(prsm_ptrs, filtered_prsm_ptrs);
  std::sort(filtered_prsm_ptrs.begin(), filtered_prsm_ptrs.end(),
            PrsmStr::cmpEValueIncProtInc);

  setProtId(filtered_prsm_ptrs);
  // find proteoform clusters and add proteoform id and intensity information
  setProteoClusterId(filtered_prsm_ptrs, prec_error_tole);
  std::sort(filtered_prsm_ptrs.begin(), filtered_prsm_ptrs.end(), 
            PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  // output
  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(filtered_prsm_ptrs);
  writer.close();
}

}

}  // namespace toppic
