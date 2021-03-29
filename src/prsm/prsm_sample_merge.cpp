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

#include "common/base/mod_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_sample_merge.hpp"

namespace toppic {

PrsmSampleMerge::PrsmSampleMerge(const std::string &db_file_name,
                                 const std::vector<std::string> &input_file_names,
                                 const std::string &output_file_name,
                                 const std::string &fix_mod,
                                 double error_tole):
    db_file_name_(db_file_name),
    input_file_names_(input_file_names),
    output_file_name_(output_file_name),
    fix_mod_(fix_mod),
    error_tole_(error_tole) {}


void PrsmSampleMerge::getPrsmClusters(PrsmStrPtrVec& prsm_ptrs,
                                      PrsmStrPtrVec2D& clusters) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    bool is_found = false;
    PrsmStrPtr cur_ptr = prsm_ptrs[i];
    for (size_t j = 0; j < clusters.size(); j++) {
      PrsmStrPtr ref_ptr = clusters[j][0];
      if (cur_ptr->getProtId() == ref_ptr->getProtId()) {
        if (std::abs(cur_ptr->getAdjustedPrecMass() - ref_ptr->getAdjustedPrecMass()) 
            <= error_tole_) {
          clusters[j].push_back(cur_ptr);
          is_found = true;
          break;
        }
      } else if (cur_ptr->getProteinMatchSeq() == ref_ptr->getProteinMatchSeq()) {
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


void PrsmSampleMerge::convertClustersToTable(PrsmStrPtrVec2D &clusters, 
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

void PrsmSampleMerge::outputTable(PrsmStrPtrVec2D &clusters,
                                  PrsmStrPtrVec2D &table_prsms,
                                  int sample_num) {
  std::ofstream file;
  file.open(output_file_name_.c_str());
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
    file << input_file_names_[i] << " abundance" << ","
        << input_file_names_[i] << " E-value" << ",";
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
        << prsm_ptr->getProteinMatchSeq() << ","
        << prsm_ptr->getVariablePtmNum() << ","
        << prsm_ptr->getUnexpectedPtmNum() << ",";
    for (int j = 0; j < sample_num; j++) {
      PrsmStrPtr sample_prsm = table_prsms[i][j];
      if (sample_prsm == nullptr) {
        file << "," << ",";
      }
      else {
        if (sample_prsm->getPrecFeatureInte() > 0) {
          file << sample_prsm->getPrecFeatureInte() << ",";
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

void PrsmSampleMerge::process() {
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name_);
  ModPtrVec fix_mod_ptr_vec = mod_util::geneFixedModList(fix_mod_);
  PrsmStrPtrVec all_prsms; 
  size_t sample_num = input_file_names_.size();
  for (size_t k = 0; k < sample_num; k++) {
    std::string input_file_name = input_file_names_[k];
    PrsmStrPtrVec prsms = PrsmReader::readAllPrsmStrsMatchSeq(input_file_name, seq_reader,
                                                              fix_mod_ptr_vec);
    for (size_t i = 0; i < prsms.size(); i++) {
      prsms[i]->setSampleId(k);
    }
    all_prsms.insert(all_prsms.end(), prsms.begin(), prsms.end());
  }
  std::sort(all_prsms.begin(), all_prsms.end(), PrsmStr::cmpEValueInc);
  PrsmStrPtrVec2D clusters; 
  getPrsmClusters(all_prsms, clusters);
  
  PrsmStrPtrVec2D table_prsms; 
  convertClustersToTable(clusters, table_prsms, sample_num);
  outputTable(clusters, table_prsms, sample_num);
}

}  // namespace toppic


