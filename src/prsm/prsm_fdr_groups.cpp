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

#include "prsm/prsm_fdr_groups.hpp"

#include <algorithm>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace toppic {

namespace prsm_fdr_groups {

// The grouping and FDR helpers are shared with prsm_fdr (prsm_fdr.hpp).

namespace {

// The FDR group of a PrSM: 0 = no mass shift and no variable PTM, 1 = no
// mass shift and variable PTMs, 2 = one mass shift, 3 = two mass shifts.
const int kGroupNum = 4;

int fdrGroup(const PrsmStrPtr& prsm_ptr) {
  int shift_num = prsm_ptr->getUnexpectedPtmNum();
  if (shift_num == 0) {
    return prsm_ptr->getVariablePtmNum() == 0 ? 0 : 1;
  }
  return shift_num == 1 ? 2 : 3;
}

}  // namespace

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
  // E-value order, which computePrsmFdr needs; the per-group target and
  // decoy lists inherit it
  std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
            PrsmStr::cmpEValueIncProtInc);

  // One FDR per search and, inside it, one per group. TopPIC's E-values come
  // from the generating function and TopMG's from MCMC sampling, so their
  // decoy distributions are not comparable and a shared count lets whichever
  // search has the smallest decoy E-value set the cutoff for both
  // (prsm_fdr::process says the same of its single FDR).
  std::map<std::string, PrsmStrPtrVec> by_tool;
  for (size_t i = 0; i < filtered_prsm_str_ptrs.size(); i++) {
    by_tool[filtered_prsm_str_ptrs[i]->getTool()].push_back(
        filtered_prsm_str_ptrs[i]);
  }

  for (auto& tool_prsms : by_tool) {
    // the PrSM level: the four groups of this search
    //   1: no mass shift, no variable PTM     3: one mass shift
    //   2: no mass shift, variable PTMs       4: two mass shifts
    PrsmStrPtrVec target_prsms[kGroupNum];
    PrsmStrPtrVec decoy_prsms[kGroupNum];
    for (size_t i = 0; i < tool_prsms.second.size(); i++) {
      PrsmStrPtr prsm_ptr = tool_prsms.second[i];
      int group = fdrGroup(prsm_ptr);
      if (prsm_fdr::isDecoy(prsm_ptr)) {
        decoy_prsms[group].push_back(prsm_ptr);
      } else {
        target_prsms[group].push_back(prsm_ptr);
      }
    }
    for (int group = 0; group < kGroupNum; group++) {
      prsm_fdr::computePrsmFdr(target_prsms[group], decoy_prsms[group]);
    }

    // The proteoform level. Targets and decoys are grouped together by
    // proteoform cluster (in E-value order, so the first PrSM of a group is
    // its best one). Each proteoform goes to the FDR group of its best PrSM
    // and is a target or a decoy according to that PrSM
    // (prsm_fdr::splitTargetDecoyProteoforms).
    PrsmStrPtrVec2D all_proteoforms =
        prsm_fdr::getProteoformGroups(tool_prsms.second);
    PrsmStrPtrVec2D group_proteoforms[kGroupNum];
    for (size_t i = 0; i < all_proteoforms.size(); i++) {
      group_proteoforms[fdrGroup(all_proteoforms[i][0])].push_back(
          all_proteoforms[i]);
    }
    for (int group = 0; group < kGroupNum; group++) {
      PrsmStrPtrVec2D target_proteoforms;
      PrsmStrPtrVec2D decoy_proteoforms;
      prsm_fdr::splitTargetDecoyProteoforms(group_proteoforms[group],
                                            target_proteoforms,
                                            decoy_proteoforms);
      prsm_fdr::computeProteoformFdr(target_proteoforms, decoy_proteoforms);
    }

    // the protein level: one target-decoy competition over this search's
    // proteoforms, as in prsm_fdr, since one protein can collect proteoforms
    // from several groups
    PrsmStrPtrVec2D target_proteoforms;
    PrsmStrPtrVec2D decoy_proteoforms;
    prsm_fdr::splitTargetDecoyProteoforms(all_proteoforms, target_proteoforms,
                                          decoy_proteoforms);
    PrsmStrPtrVec2D target_proteins =
        prsm_fdr::getProteinGroups(target_proteoforms);
    PrsmStrPtrVec2D decoy_proteins =
        prsm_fdr::getProteinGroups(decoy_proteoforms);
    prsm_fdr::computeProteinFdr(target_proteins, decoy_proteins);
  }

  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);

  // every PrSM, targets and decoys, in spectrum order
  std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
            PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  writer.writeVector(filtered_prsm_str_ptrs);
  writer.close();
}

}  // namespace prsm_fdr_groups

}  // namespace toppic
