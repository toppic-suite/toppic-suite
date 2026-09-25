// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "prsm/prsm_simple_cluster.hpp"

#include <algorithm>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace toppic {

namespace prsm_simple_cluster {

PrsmStrPtrVec2D setProtId(PrsmStrPtrVec& prsm_ptrs) {
  PrsmStrPtrVec2D proteins;
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
  return proteins;
}

// Two proteoforms of one protein belong to the same cluster only if their
// sequences overlap by at least this fraction of the shorter one; without the
// check, two PrSMs of the same protein and mass covering disjoint regions of
// the protein (an N-terminal and a C-terminal fragment of equal mass) are
// clustered together.
const double kMinOverlapRatio = 0.7;

// Start and end positions are inclusive, so a one-residue proteoform has
// length 1.
bool isOverlapMatch(const PrsmStrPtr& a, const PrsmStrPtr& b) {
  int a_start = a->getProteoformStartPos();
  int a_end = a->getProteoformEndPos();
  int b_start = b->getProteoformStartPos();
  int b_end = b->getProteoformEndPos();
  int overlap = std::min(a_end, b_end) - std::max(a_start, b_start) + 1;
  if (overlap <= 0) {
    return false;
  }
  int shorter_len = std::min(a_end - a_start + 1, b_end - b_start + 1);
  if (shorter_len <= 0) {
    return false;
  }
  return overlap >= kMinOverlapRatio * shorter_len;
}

void setClusterId(PrsmStrPtrVec2D& proteins, bool is_ppm_error,
                  double error_tole) {
  PrsmStrPtrVec2D clusters;
  double dalton_error_tole = error_tole;
  for (size_t i = 0; i < proteins.size(); i++) {
    PrsmStrPtrVec2D protein_clusters;
    for (size_t j = 0; j < proteins[i].size(); j++) {
      bool is_found = false;
      PrsmStrPtr cur_prsm = proteins[i][j];
      if (is_ppm_error) {
        dalton_error_tole = cur_prsm->getOriPrecMass() * error_tole;
      }
      for (size_t m = 0; m < protein_clusters.size(); m++) {
        PrsmStrPtr ref_prsm = protein_clusters[m][0];
        if (PrsmStr::isSimpleMatch(cur_prsm, ref_prsm, dalton_error_tole) &&
            isOverlapMatch(cur_prsm, ref_prsm)) {
          protein_clusters[m].push_back(cur_prsm);
          is_found = true;
          break;
        }
      }
      if (!is_found) {
        PrsmStrPtrVec new_clusters;
        new_clusters.push_back(cur_prsm);
        protein_clusters.push_back(new_clusters);
      }
    }
    clusters.insert(std::end(clusters), std::begin(protein_clusters),
                    std::end(protein_clusters));
  }

  for (size_t i = 0; i < clusters.size(); i++) {
    for (size_t j = 0; j < clusters[i].size(); j++) {
      clusters[i][j]->setProteoClusterId(i);
    }
  }
}

void process(const std::string& db_file_name, const std::string& spec_file_name,
             const std::string& input_file_ext,
             const ModPtrVec& fix_mod_ptr_vec,
             const std::string& output_file_ext, bool is_ppm_error,
             double error_tole) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;
  LOG_DEBUG("Reading prsm strings started");
  PrsmStrPtrVec prsm_ptrs = prsm_reader_util::readAllPrsmStrs(input_file_name);
  LOG_DEBUG("Reading prsm strings finished");
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpEValueIncProtInc);
  PrsmStrPtrVec2D protein_prsms = setProtId(prsm_ptrs);
  setClusterId(protein_prsms, is_ppm_error, error_tole);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(),
       PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  // output
  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prsm_simple_cluster

}  // namespace toppic
