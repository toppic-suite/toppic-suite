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

#include "common/base/mod_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_sample_merge.hpp"

namespace toppic {

namespace prsm_sample_merge {

void getPrsmClusters(PrsmStrPtrVec& prsm_ptrs,
                     PrsmStrPtrVec2D& clusters, 
                     double error_tole) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    bool is_found = false;
    PrsmStrPtr cur_ptr = prsm_ptrs[i];
    for (size_t j = 0; j < clusters.size(); j++) {
      PrsmStrPtr ref_ptr = clusters[j][0];
      if (cur_ptr->getProtId() == ref_ptr->getProtId()) {
        if (std::abs(cur_ptr->getAdjustedPrecMass() - ref_ptr->getAdjustedPrecMass()) 
            <= error_tole) {
          clusters[j].push_back(cur_ptr);
          is_found = true;
          break;
        }
      } else if (cur_ptr->getProteoformMatchSeq() == ref_ptr->getProteoformMatchSeq()) {
        clusters[j].push_back(cur_ptr);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmStrPtrVec new_clusters;
      new_clusters.push_back(prsm_ptrs[i]);
      clusters.push_back(new_clusters);
    }
  }
}

void convertClustersToTable(PrsmStrPtrVec2D &clusters, 
                            PrsmStrPtrVec2D &table_prsms,
                            int sample_num) {
  int cluster_num = clusters.size();
  PrsmStrPtr null_prsm = nullptr;
  for (int i = 0; i < cluster_num; i++) {
    PrsmStrPtrVec cur_prsms(sample_num);
    for (size_t j = 0; j < clusters[i].size(); j++) {
      PrsmStrPtr cur_prsm = clusters[i][j];
      int sample_id = cur_prsm->getSampleId();
      cur_prsms[sample_id] = cur_prsm;
    }
    table_prsms.push_back(cur_prsms);
  }
}

void outputTable(PrsmStrPtrVec2D &clusters,
                 PrsmStrPtrVec2D &table_prsms,
                 int sample_num,
                 const std::vector<std::string> &input_file_names,
                 const std::string &output_file_name) {
  std::ofstream file;
  file.open(output_file_name.c_str());
  // write title
  file << "Protein accession" << ","
      << "Protein description" << ","
      << "Adjusted precursor mass" << ","
      << "First residue" << ","
      << "Last residue" << ","
      << "Proteoform" << ","
      << "#variable PTMs" << ","
      << "#unexpected modifications" << ",";

  for (int i = 0; i < sample_num; i++) {
    file << input_file_names[i] << " abundance" << ","
        << input_file_names[i] << " E-value" << ",";
  }
  file << std::endl;
  int cluster_num = clusters.size();
  for (int i = 0; i < cluster_num; i++) {
    PrsmStrPtr prsm_ptr = clusters[i][0];
    file << prsm_ptr->getSeqName() << ","
        << "\"" << prsm_ptr->getSeqDesc() << "\"" << ","
        << prsm_ptr->getAdjustedPrecMass() << ","
        << (prsm_ptr->getProteoformStartPos() + 1) << ","
        << (prsm_ptr->getProteoformEndPos() + 1) << ","
        << prsm_ptr->getProteoformMatchSeq() << ","
        << prsm_ptr->getVariablePtmNum() << ","
        << prsm_ptr->getUnexpectedPtmNum() << ",";
    for (int j = 0; j < sample_num; j++) {
      PrsmStrPtr sample_prsm = table_prsms[i][j];
      if (sample_prsm == nullptr) {
        file << "," << ",";
      }
      else {
        if (sample_prsm->getFracFeatureInte() > 0) {
          file << sample_prsm->getFracFeatureInte() << ",";
        } else {
          file << "-" << ",";
        }
        file << sample_prsm->getEValue() << ",";
      }
    }
    file << std::endl;
  }
  file.close();
}

void process(const std::vector<std::string> &input_file_names,
             const std::string &output_file_name,
             double error_tole) {
  PrsmStrPtrVec all_prsms; 
  size_t sample_num = input_file_names.size();
  for (size_t k = 0; k < sample_num; k++) {
    std::string input_file_name = input_file_names[k];
    PrsmStrPtrVec prsms = prsm_reader_util::readAllPrsmStrsMatchSeq(input_file_name);
    /*
    for (size_t i = 0; i < prsms.size(); i++) {
      prsms[i]->setSampleId(k);
    }
    */
    all_prsms.insert(all_prsms.end(), prsms.begin(), prsms.end());
  }
  std::sort(all_prsms.begin(), all_prsms.end(), PrsmStr::cmpEValueIncProtInc);
  PrsmStrPtrVec2D clusters; 
  getPrsmClusters(all_prsms, clusters, error_tole);
  
  PrsmStrPtrVec2D table_prsms; 
  convertClustersToTable(clusters, table_prsms, sample_num);
  outputTable(clusters, table_prsms, sample_num, input_file_names, output_file_name);
}

}

}  // namespace toppic


