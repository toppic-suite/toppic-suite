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

#include "prsm/prsm_fdr.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace toppic {

namespace prsm_fdr {

// Group the PrSMs by proteoform cluster (see prsm_feature_cluster /
// prsm_simple_cluster). The input is sorted by E-value, so the first PrSM of
// a group is the proteoform's best one (lowest E-value) and the groups come
// out ordered by that E-value.
PrsmStrPtrVec2D getGroups(const PrsmStrPtrVec& prsm_ptrs) {
  PrsmStrPtrVec2D results;
  std::map<int, size_t> group_idx;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int cluster_id = prsm_ptrs[i]->getProteoClusterId();
    auto it = group_idx.find(cluster_id);
    if (it == group_idx.end()) {
      group_idx[cluster_id] = results.size();
      results.push_back(PrsmStrPtrVec{prsm_ptrs[i]});
    } else {
      results[it->second].push_back(prsm_ptrs[i]);
    }
  }
  return results;
}

// Group the proteoform groups of getGroups by protein cluster (see
// prsm_prot_cluster), concatenating the PrSMs of the proteoforms of one
// protein. The input groups are ordered by their best E-value, so the first
// PrSM of a protein group is the protein's best proteoform's best PrSM
// (lowest E-value) and the protein groups come out ordered by that E-value.
PrsmStrPtrVec2D getProteinGroups(const PrsmStrPtrVec2D& proteoforms) {
  PrsmStrPtrVec2D results;
  std::map<int, size_t> group_idx;
  for (size_t i = 0; i < proteoforms.size(); i++) {
    int cluster_id = proteoforms[i][0]->getProtClusterId();
    auto it = group_idx.find(cluster_id);
    if (it == group_idx.end()) {
      group_idx[cluster_id] = results.size();
      results.push_back(proteoforms[i]);
    } else {
      results[it->second].insert(results[it->second].end(),
                                 proteoforms[i].begin(), proteoforms[i].end());
    }
  }
  return results;
}

// Split the proteoform groups of getGroups (targets and decoys grouped
// together, ordered by their best E-value) by their best PrSM: a proteoform
// is a target proteoform if its best (lowest E-value) PrSM is a target PrSM,
// a decoy proteoform otherwise. Both outputs keep the best-E-value order, as
// computeProteoformFdr needs.
void splitTargetDecoyProteoforms(const PrsmStrPtrVec2D& all_proteoforms,
                                 PrsmStrPtrVec2D& target_proteoforms,
                                 PrsmStrPtrVec2D& decoy_proteoforms) {
  for (size_t i = 0; i < all_proteoforms.size(); i++) {
    if (isDecoy(all_proteoforms[i][0])) {
      decoy_proteoforms.push_back(all_proteoforms[i]);
    } else {
      target_proteoforms.push_back(all_proteoforms[i]);
    }
  }
}

namespace {

double fdrRatio(int n_decoy, int n_target) {
  double fdr = static_cast<double>(n_decoy) / static_cast<double>(n_target);
  if (fdr > 1) {
    fdr = 1.0;
  }
  return fdr;
}

}  // namespace

void computeFdr(PrsmStrPtrVec& target_ptrs, PrsmStrPtrVec& decoy_ptrs) {
  int n_decoy = 0;
  for (size_t i = 0; i < target_ptrs.size(); i++) {
    int n_target = i + 1;
    double target_evalue = target_ptrs[i]->getEValue();
    for (size_t j = n_decoy; j < decoy_ptrs.size(); j++) {
      if (decoy_ptrs[j]->getEValue() <= target_evalue) {
        n_decoy++;
        double fdr = fdrRatio(n_decoy, n_target);
        decoy_ptrs[j]->setFdr(fdr);
      } else {
        break;
      }
    }
    double fdr = fdrRatio(n_decoy, n_target);
    target_ptrs[i]->setFdr(fdr);
  }
}

void computeProteoformFdr(PrsmStrPtrVec2D& target_proteoforms,
                          PrsmStrPtrVec2D& decoy_proteoforms) {
  int n_decoy = 0;
  for (size_t i = 0; i < target_proteoforms.size(); i++) {
    int n_target = i + 1;
    double target_evalue = target_proteoforms[i][0]->getEValue();
    for (size_t j = n_decoy; j < decoy_proteoforms.size(); j++) {
      if (decoy_proteoforms[j][0]->getEValue() <= target_evalue) {
        n_decoy++;
        double fdr = fdrRatio(n_decoy, n_target);
        for (size_t k = 0; k < decoy_proteoforms[j].size(); k++) {
          decoy_proteoforms[j][k]->setProteoformFdr(fdr);
        }
      } else {
        break;
      }
    }
    double fdr = fdrRatio(n_decoy, n_target);
    for (size_t k = 0; k < target_proteoforms[i].size(); k++) {
      target_proteoforms[i][k]->setProteoformFdr(fdr);
    }
  }
}

// Protein-level FDR: each protein is represented by its best proteoform (the
// PrSM with the lowest E-value, the first of its group); target and decoy
// proteins are walked in E-value order and the FDR at a target protein is
// (#decoy proteins with a better or equal E-value) / (#target proteins so
// far). The value is assigned to every PrSM of the protein.
void computeProteinFdr(PrsmStrPtrVec2D& target_proteins,
                       PrsmStrPtrVec2D& decoy_proteins) {
  int n_decoy = 0;
  for (size_t i = 0; i < target_proteins.size(); i++) {
    int n_target = i + 1;
    double target_evalue = target_proteins[i][0]->getEValue();
    for (size_t j = n_decoy; j < decoy_proteins.size(); j++) {
      if (decoy_proteins[j][0]->getEValue() <= target_evalue) {
        n_decoy++;
        double fdr = fdrRatio(n_decoy, n_target);
        for (size_t k = 0; k < decoy_proteins[j].size(); k++) {
          decoy_proteins[j][k]->setProteinFdr(fdr);
        }
      } else {
        break;
      }
    }
    double fdr = fdrRatio(n_decoy, n_target);
    for (size_t k = 0; k < target_proteins[i].size(); k++) {
      target_proteins[i][k]->setProteinFdr(fdr);
    }
  }
  // decoy proteins ranked below every target protein
  for (size_t j = n_decoy; j < decoy_proteins.size(); j++) {
    for (size_t k = 0; k < decoy_proteins[j].size(); k++) {
      decoy_proteins[j][k]->setProteinFdr(1.0);
    }
  }
}

void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;

  // targets and decoys together, without the zero-E-value PrSMs
  PrsmStrPtrVec prsm_str_ptrs =
      prsm_reader_util::readAllPrsmStrs(input_file_name);
  PrsmStrPtrVec filtered_prsm_str_ptrs;
  for (size_t i = 0; i < prsm_str_ptrs.size(); i++) {
    if (prsm_str_ptrs[i]->getEValue() == 0.0) {
      LOG_ERROR("toppic::PRSMFdr zero E value is reported");
    } else {
      filtered_prsm_str_ptrs.push_back(prsm_str_ptrs[i]);
    }
  }
  // E-value order, which computeFdr and getGroups need; the target and decoy
  // lists inherit it
  std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
            PrsmStr::cmpEValueIncProtInc);
  PrsmStrPtrVec target_ptrs;
  PrsmStrPtrVec decoy_ptrs;
  for (size_t i = 0; i < filtered_prsm_str_ptrs.size(); i++) {
    if (isDecoy(filtered_prsm_str_ptrs[i])) {
      decoy_ptrs.push_back(filtered_prsm_str_ptrs[i]);
    } else {
      target_ptrs.push_back(filtered_prsm_str_ptrs[i]);
    }
  }
  computeFdr(target_ptrs, decoy_ptrs);

  PrsmStrPtrVec2D all_proteoforms = getGroups(filtered_prsm_str_ptrs);
  PrsmStrPtrVec2D target_proteoforms;
  PrsmStrPtrVec2D decoy_proteoforms;
  splitTargetDecoyProteoforms(all_proteoforms, target_proteoforms,
                              decoy_proteoforms);
  computeProteoformFdr(target_proteoforms, decoy_proteoforms);

  PrsmStrPtrVec2D target_proteins = getProteinGroups(target_proteoforms);
  PrsmStrPtrVec2D decoy_proteins = getProteinGroups(decoy_proteoforms);
  computeProteinFdr(target_proteins, decoy_proteins);

  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
            PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  std::cout << "target_ptrs: " << target_ptrs.size()
            << ", decoy_ptrs: " << decoy_ptrs.size() << std::endl;
  writer.writeVector(filtered_prsm_str_ptrs);
  writer.close();
}

}  // namespace prsm_fdr

}  // namespace toppic
