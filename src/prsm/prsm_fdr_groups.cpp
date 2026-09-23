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

  // Group 1: No mass shift and 0 variable PTMS
  PrsmStrPtrVec target_prsms1;
  PrsmStrPtrVec decoy_prsms1;
  // Group2: No mass shift and >0 variable PTMs
  PrsmStrPtrVec target_prsms2;
  PrsmStrPtrVec decoy_prsms2;
  // Group3: 1 mass shift
  PrsmStrPtrVec target_prsms3;
  PrsmStrPtrVec decoy_prsms3;
  // Group4: 2 mass shift
  PrsmStrPtrVec target_prsms4;
  PrsmStrPtrVec decoy_prsms4;
  for (size_t i = 0; i < filtered_prsm_str_ptrs.size(); i++) {
    if (prsm_fdr::isDecoy(filtered_prsm_str_ptrs[i])) {
      // 0 mass shift
      if (filtered_prsm_str_ptrs[i]->getUnexpectedPtmNum() == 0) {
        // 0 variable PTMs
        if (filtered_prsm_str_ptrs[i]->getVariablePtmNum() == 0) {
          decoy_prsms1.push_back(filtered_prsm_str_ptrs[i]);
        }
        //>0 variable PTMs
        else {
          decoy_prsms2.push_back(filtered_prsm_str_ptrs[i]);
        }
      }
      // 1 mass shift
      else if (filtered_prsm_str_ptrs[i]->getUnexpectedPtmNum() == 1) {
        decoy_prsms3.push_back(filtered_prsm_str_ptrs[i]);
      }
      // 2 mass shift
      else {
        decoy_prsms4.push_back(filtered_prsm_str_ptrs[i]);
      }
    } else {
      // 0 mass shift
      if (filtered_prsm_str_ptrs[i]->getUnexpectedPtmNum() == 0) {
        // 0 variable PTMs
        if (filtered_prsm_str_ptrs[i]->getVariablePtmNum() == 0) {
          target_prsms1.push_back(filtered_prsm_str_ptrs[i]);
        }
        //>0 variable PTMs
        else {
          target_prsms2.push_back(filtered_prsm_str_ptrs[i]);
        }
      }
      // 1 mass shift
      else if (filtered_prsm_str_ptrs[i]->getUnexpectedPtmNum() == 1) {
        target_prsms3.push_back(filtered_prsm_str_ptrs[i]);
      }
      // 2 mass shift
      else {
        target_prsms4.push_back(filtered_prsm_str_ptrs[i]);
      }
    }
  }

  prsm_fdr::computePrsmFdr(target_prsms1, decoy_prsms1);
  prsm_fdr::computePrsmFdr(target_prsms2, decoy_prsms2);
  prsm_fdr::computePrsmFdr(target_prsms3, decoy_prsms3);
  prsm_fdr::computePrsmFdr(target_prsms4, decoy_prsms4);

  // Targets and decoys are grouped together by proteoform cluster (in
  // E-value order, so the first PrSM of a group is its best one). Each
  // proteoform goes to the FDR group of its best PrSM and is a target or a
  // decoy according to that PrSM (prsm_fdr::splitTargetDecoyProteoforms).
  PrsmStrPtrVec2D all_proteoforms =
      prsm_fdr::getProteoformGroups(filtered_prsm_str_ptrs);
  PrsmStrPtrVec2D group_proteoforms[4];
  for (size_t i = 0; i < all_proteoforms.size(); i++) {
    group_proteoforms[fdrGroup(all_proteoforms[i][0])].push_back(
        all_proteoforms[i]);
  }
  PrsmStrPtrVec2D target1_proteoforms;
  PrsmStrPtrVec2D decoy1_proteoforms;
  PrsmStrPtrVec2D target2_proteoforms;
  PrsmStrPtrVec2D decoy2_proteoforms;
  PrsmStrPtrVec2D target3_proteoforms;
  PrsmStrPtrVec2D decoy3_proteoforms;
  PrsmStrPtrVec2D target4_proteoforms;
  PrsmStrPtrVec2D decoy4_proteoforms;
  prsm_fdr::splitTargetDecoyProteoforms(group_proteoforms[0],
                                        target1_proteoforms,
                                        decoy1_proteoforms);
  prsm_fdr::splitTargetDecoyProteoforms(group_proteoforms[1],
                                        target2_proteoforms,
                                        decoy2_proteoforms);
  prsm_fdr::splitTargetDecoyProteoforms(group_proteoforms[2],
                                        target3_proteoforms,
                                        decoy3_proteoforms);
  prsm_fdr::splitTargetDecoyProteoforms(group_proteoforms[3],
                                        target4_proteoforms,
                                        decoy4_proteoforms);

  prsm_fdr::computeProteoformFdr(target1_proteoforms, decoy1_proteoforms);
  prsm_fdr::computeProteoformFdr(target2_proteoforms, decoy2_proteoforms);
  prsm_fdr::computeProteoformFdr(target3_proteoforms, decoy3_proteoforms);
  prsm_fdr::computeProteoformFdr(target4_proteoforms, decoy4_proteoforms);

  // protein-level FDR over all proteoforms, as in prsm_fdr: the per-group
  // FDRs above leave the protein level to a single target-decoy competition,
  // since one protein can collect proteoforms from several groups
  PrsmStrPtrVec2D target_proteoforms;
  PrsmStrPtrVec2D decoy_proteoforms;
  prsm_fdr::splitTargetDecoyProteoforms(all_proteoforms, target_proteoforms,
                                        decoy_proteoforms);
  PrsmStrPtrVec2D target_proteins =
      prsm_fdr::getProteinGroups(target_proteoforms);
  PrsmStrPtrVec2D decoy_proteins =
      prsm_fdr::getProteinGroups(decoy_proteoforms);
  prsm_fdr::computeProteinFdr(target_proteins, decoy_proteins);

  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);

  // every PrSM, targets and decoys, in spectrum order
  std::sort(filtered_prsm_str_ptrs.begin(), filtered_prsm_str_ptrs.end(),
            PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  writer.writeVector(filtered_prsm_str_ptrs);
  std::cout << "Group 1 target_prsms: " << target_prsms1.size()
            << ", Group 1 decoy_prsms: " << decoy_prsms1.size() << std::endl;
  std::cout << "Group 2 target_prsms: " << target_prsms2.size()
            << ", Group 2 decoy_prsms: " << decoy_prsms2.size() << std::endl;
  std::cout << "Group 3 target_prsms: " << target_prsms3.size()
            << ", Group 3 decoy_prsms: " << decoy_prsms3.size() << std::endl;
  std::cout << "Group 4 target_prsms: " << target_prsms4.size()
            << ", Group 4 decoy_prsms: " << decoy_prsms4.size() << std::endl;
  writer.close();
}

}  // namespace prsm_fdr_groups

}  // namespace toppic
