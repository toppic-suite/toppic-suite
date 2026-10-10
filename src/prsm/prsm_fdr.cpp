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
#include <iostream>
#include <map>
#include <memory>
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
PrsmStrPtrVec2D getProteoformGroups(const PrsmStrPtrVec& prsm_ptrs) {
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

// Group the proteoform groups of getProteoformGroups by protein cluster (see
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

// Split the proteoform groups of getProteoformGroups (targets and decoys
// grouped together, ordered by their best E-value) by their best PrSM: a
// proteoform is a target proteoform if its best (lowest E-value) PrSM is a
// target PrSM, a decoy proteoform otherwise. Both outputs keep the
// best-E-value order, as computeProteoformFdr needs.
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

void computePrsmFdr(PrsmStrPtrVec& target_prsms, PrsmStrPtrVec& decoy_prsms,
                    FdrKey key) {
  int n_decoy = 0;
  for (size_t i = 0; i < target_prsms.size(); i++) {
    int n_target = i + 1;
    double target_evalue = key(target_prsms[i]);
    for (size_t j = n_decoy; j < decoy_prsms.size(); j++) {
      if (key(decoy_prsms[j]) <= target_evalue) {
        n_decoy++;
        double fdr = fdrRatio(n_decoy, n_target);
        decoy_prsms[j]->setFdr(fdr);
      } else {
        break;
      }
    }
    double fdr = fdrRatio(n_decoy, n_target);
    target_prsms[i]->setFdr(fdr);
  }
}

void computeProteoformFdr(PrsmStrPtrVec2D& target_proteoforms,
                          PrsmStrPtrVec2D& decoy_proteoforms, FdrKey key) {
  int n_decoy = 0;
  for (size_t i = 0; i < target_proteoforms.size(); i++) {
    int n_target = i + 1;
    double target_evalue = key(target_proteoforms[i][0]);
    for (size_t j = n_decoy; j < decoy_proteoforms.size(); j++) {
      if (key(decoy_proteoforms[j][0]) <= target_evalue) {
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
  // A decoy proteoform (a cluster whose best PrSM is a decoy) may hold
  // target PrSMs of other spectra of the cluster. They are not a target
  // proteoform and must not inherit the decoy proteoform's FDR, which would
  // let them pass the proteoform-level cutoff (toppic_test experiment 40):
  // they get FDR 1. (The decoy PrSMs of a target proteoform keep the target
  // proteoform's FDR; they are only reported with -K.)
  for (size_t j = 0; j < decoy_proteoforms.size(); j++) {
    for (size_t k = 0; k < decoy_proteoforms[j].size(); k++) {
      if (!isDecoy(decoy_proteoforms[j][k])) {
        decoy_proteoforms[j][k]->setProteoformFdr(1.0);
      }
    }
  }
}

// Protein-level FDR: each protein is represented by its best proteoform (the
// PrSM with the lowest E-value, the first of its group); target and decoy
// proteins are walked in E-value order and the FDR at a target protein is
// (#decoy proteins with a better or equal E-value) / (#target proteins so
// far). The value is assigned to every PrSM of the protein.
void computeProteinFdr(PrsmStrPtrVec2D& target_proteins,
                       PrsmStrPtrVec2D& decoy_proteins, FdrKey key) {
  int n_decoy = 0;
  for (size_t i = 0; i < target_proteins.size(); i++) {
    int n_target = i + 1;
    double target_evalue = key(target_proteins[i][0]);
    for (size_t j = n_decoy; j < decoy_proteins.size(); j++) {
      if (key(decoy_proteins[j][0]) <= target_evalue) {
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
  // the target PrSMs of a decoy protein (see computeProteoformFdr) are not
  // a target protein: FDR 1
  for (size_t j = 0; j < decoy_proteins.size(); j++) {
    for (size_t k = 0; k < decoy_proteins[j].size(); k++) {
      if (!isDecoy(decoy_proteins[j][k])) {
        decoy_proteins[j][k]->setProteinFdr(1.0);
      }
    }
  }
}

namespace {

// the targets within a level's cutoff, and the search's target and decoy
// PrSMs
struct PassCounts {
  int prsm = 0;
  int proteoform = 0;
  int protein = 0;
  size_t target_prsm_num = 0;
  size_t decoy_prsm_num = 0;
};

// The spectrum-, proteoform- and protein-level FDRs of one search's PrSMs,
// which are sorted by the key (the order getProteoformGroups needs), and the
// targets within the cutoffs: target PrSMs, and target proteoforms and
// proteins (groups of a target best PrSM) by their FDR.
PassCounts computeFdrs(const PrsmStrPtrVec& prsms, FdrKey key,
                       const FdrCutoffs& cutoffs) {
  PrsmStrPtrVec target_prsms;
  PrsmStrPtrVec decoy_prsms;
  for (size_t i = 0; i < prsms.size(); i++) {
    if (isDecoy(prsms[i])) {
      decoy_prsms.push_back(prsms[i]);
    } else {
      target_prsms.push_back(prsms[i]);
    }
  }
  computePrsmFdr(target_prsms, decoy_prsms, key);

  PrsmStrPtrVec2D all_proteoforms = getProteoformGroups(prsms);
  PrsmStrPtrVec2D target_proteoforms;
  PrsmStrPtrVec2D decoy_proteoforms;
  splitTargetDecoyProteoforms(all_proteoforms, target_proteoforms,
                              decoy_proteoforms);
  computeProteoformFdr(target_proteoforms, decoy_proteoforms, key);

  PrsmStrPtrVec2D target_proteins = getProteinGroups(target_proteoforms);
  PrsmStrPtrVec2D decoy_proteins = getProteinGroups(decoy_proteoforms);
  computeProteinFdr(target_proteins, decoy_proteins, key);

  PassCounts counts;
  counts.target_prsm_num = target_prsms.size();
  counts.decoy_prsm_num = decoy_prsms.size();
  for (const PrsmStrPtr& prsm_ptr : target_prsms) {
    if (prsm_ptr->getFdr() <= cutoffs.prsm) {
      counts.prsm++;
    }
  }
  for (const PrsmStrPtrVec& group : target_proteoforms) {
    if (group[0]->getProteoformFdr() <= cutoffs.proteoform) {
      counts.proteoform++;
    }
  }
  for (const PrsmStrPtrVec& group : target_proteins) {
    if (group[0]->getProteinFdr() <= cutoffs.protein) {
      counts.protein++;
    }
  }
  return counts;
}

}  // namespace

void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext, bool use_rescore,
             const FdrCutoffs& cutoffs, bool rescore_fallback) {
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

  // The E-value round of a rescored FDR computation runs on copies of the
  // PrSMs, so that the rescore round below writes the PrSMs as a single
  // rescored round does; the copies are sorted from the input order as the
  // E-value computation sorts them (the same order without the rescoring).
  bool two_rounds = use_rescore && rescore_fallback;
  std::map<std::string, PrsmStrPtrVec> evalue_by_tool;
  std::map<const PrsmStr*, PrsmStrPtr> evalue_copies;
  if (two_rounds) {
    PrsmStrPtrVec copies;
    for (const PrsmStrPtr& prsm_ptr : filtered_prsm_str_ptrs) {
      PrsmStrPtr copy = std::make_shared<PrsmStr>(*prsm_ptr);
      copies.push_back(copy);
      evalue_copies[prsm_ptr.get()] = copy;
    }
    std::sort(copies.begin(), copies.end(), PrsmStr::cmpEValueIncProtInc);
    for (const PrsmStrPtr& copy : copies) {
      evalue_by_tool[copy->getTool()].push_back(copy);
    }
  }

  // E-value order (with the rescoring the rescoring step's score order),
  // which computePrsmFdr and getProteoformGroups need; the target and decoy
  // lists inherit it
  FdrKey key = eValueKey;
  if (use_rescore) {
    key = rescoreKey;
    std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
              PrsmStr::cmpRescoreDecProtInc);
  } else {
    std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
              PrsmStr::cmpEValueIncProtInc);
  }

  // One FDR per search. TopPIC's E-values come from the generating function
  // and TopMG's from MCMC sampling, so their decoy distributions are not
  // comparable: mixed in one target and decoy count, whichever search
  // produces the smallest decoy E-value sets the cutoff for both. Each
  // search's PrSMs are therefore counted against its own decoys, and the
  // cutoffs of TopPIC_report then keep the PrSMs that pass in either - the
  // passing sets merged.
  std::map<std::string, PrsmStrPtrVec> by_tool;
  for (size_t i = 0; i < filtered_prsm_str_ptrs.size(); i++) {
    by_tool[filtered_prsm_str_ptrs[i]->getTool()].push_back(
        filtered_prsm_str_ptrs[i]);
  }

  for (auto& tool_prsms : by_tool) {
    PassCounts evalue_counts;
    if (two_rounds) {
      evalue_counts =
          computeFdrs(evalue_by_tool[tool_prsms.first], eValueKey, cutoffs);
    }
    PassCounts counts = computeFdrs(tool_prsms.second, key, cutoffs);

    std::cout << tool_prsms.first << " target_prsms: "
              << counts.target_prsm_num
              << ", decoy_prsms: " << counts.decoy_prsm_num << std::endl;

    if (two_rounds) {
      // at each level the E-value's FDRs when they give more targets
      bool evalue_prsm = evalue_counts.prsm > counts.prsm;
      bool evalue_form = evalue_counts.proteoform > counts.proteoform;
      bool evalue_prot = evalue_counts.protein > counts.protein;
      for (const PrsmStrPtr& prsm_ptr : tool_prsms.second) {
        const PrsmStrPtr& copy = evalue_copies[prsm_ptr.get()];
        if (evalue_prsm && copy->getFdr() != prsm_ptr->getFdr()) {
          prsm_ptr->setFdr(copy->getFdr());
        }
        if (evalue_form &&
            copy->getProteoformFdr() != prsm_ptr->getProteoformFdr()) {
          prsm_ptr->setProteoformFdr(copy->getProteoformFdr());
        }
        if (evalue_prot &&
            copy->getProteinFdr() != prsm_ptr->getProteinFdr()) {
          prsm_ptr->setProteinFdr(copy->getProteinFdr());
        }
      }
      auto by = [](bool evalue) { return evalue ? "E-value" : "rescore"; };
      std::cout << "FDR computation - " << tool_prsms.first
                << ": targets within the cutoffs by E-value / rescore: PrSMs "
                << evalue_counts.prsm << " / " << counts.prsm
                << ", proteoforms " << evalue_counts.proteoform << " / "
                << counts.proteoform << ", proteins " << evalue_counts.protein
                << " / " << counts.protein << "; FDRs by " << by(evalue_prsm)
                << ", " << by(evalue_form) << ", " << by(evalue_prot) << "."
                << std::endl;
    }
  }

  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
            PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  writer.writeVector(filtered_prsm_str_ptrs);
  writer.close();
}

}  // namespace prsm_fdr

}  // namespace toppic
