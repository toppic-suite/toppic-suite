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

void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;

  PrsmStrPtrVec prsm_str_ptrs =
      prsm_reader_util::readAllPrsmStrs(input_file_name);
  PrsmStrPtrVec all_prsms1;
  PrsmStrPtrVec all_prsms2;
  PrsmStrPtrVec all_prsms3;
  PrsmStrPtrVec all_prsms4;

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
  for (size_t i = 0; i < prsm_str_ptrs.size(); i++) {
    if (prsm_str_ptrs[i]->getEValue() == 0.0) {
      LOG_ERROR("toppic::PRSMFdr zero E value is reported");
    } else {
      if (prsm_fdr::isDecoy(prsm_str_ptrs[i])) {
        // 0 mass shift
        if (prsm_str_ptrs[i]->getUnexpectedPtmNum() == 0) {
          // 0 variable PTMs
          if (prsm_str_ptrs[i]->getVariablePtmNum() == 0) {
            decoy_prsms1.push_back(prsm_str_ptrs[i]);
          }
          //>0 variable PTMs
          else {
            decoy_prsms2.push_back(prsm_str_ptrs[i]);
          }
        }
        // 1 mass shift
        else if (prsm_str_ptrs[i]->getUnexpectedPtmNum() == 1) {
          decoy_prsms3.push_back(prsm_str_ptrs[i]);
        }
        // 2 mass shift
        else {
          decoy_prsms4.push_back(prsm_str_ptrs[i]);
        }
      } else {
        // 0 mass shift
        if (prsm_str_ptrs[i]->getUnexpectedPtmNum() == 0) {
          // 0 variable PTMs
          if (prsm_str_ptrs[i]->getVariablePtmNum() == 0) {
            target_prsms1.push_back(prsm_str_ptrs[i]);
          }
          //>0 variable PTMs
          else {
            target_prsms2.push_back(prsm_str_ptrs[i]);
          }
        }
        // 1 mass shift
        else if (prsm_str_ptrs[i]->getUnexpectedPtmNum() == 1) {
          target_prsms3.push_back(prsm_str_ptrs[i]);
        }
        // 2 mass shift
        else {
          target_prsms4.push_back(prsm_str_ptrs[i]);
        }
      }
    }
  }
  std::sort(target_prsms1.begin(), target_prsms1.end(),
            PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_prsms1.begin(), decoy_prsms1.end(),
            PrsmStr::cmpEValueIncProtInc);
  std::sort(target_prsms2.begin(), target_prsms2.end(),
            PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_prsms2.begin(), decoy_prsms2.end(),
            PrsmStr::cmpEValueIncProtInc);
  std::sort(target_prsms3.begin(), target_prsms3.end(),
            PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_prsms3.begin(), decoy_prsms3.end(),
            PrsmStr::cmpEValueIncProtInc);
  std::sort(target_prsms4.begin(), target_prsms4.end(),
            PrsmStr::cmpEValueIncProtInc);
  std::sort(decoy_prsms4.begin(), decoy_prsms4.end(),
            PrsmStr::cmpEValueIncProtInc);

  prsm_fdr::computePrsmFdr(target_prsms1, decoy_prsms1);
  prsm_fdr::computePrsmFdr(target_prsms2, decoy_prsms2);
  prsm_fdr::computePrsmFdr(target_prsms3, decoy_prsms3);
  prsm_fdr::computePrsmFdr(target_prsms4, decoy_prsms4);

  // Within each group, targets and decoys are grouped together by proteoform
  // cluster in E-value order, and a proteoform is a target or a decoy
  // according to its best PrSM (prsm_fdr::splitTargetDecoyProteoforms).
  PrsmStrPtrVec2D target1_proteoforms;
  PrsmStrPtrVec2D decoy1_proteoforms;
  PrsmStrPtrVec2D target2_proteoforms;
  PrsmStrPtrVec2D decoy2_proteoforms;
  PrsmStrPtrVec2D target3_proteoforms;
  PrsmStrPtrVec2D decoy3_proteoforms;
  PrsmStrPtrVec2D target4_proteoforms;
  PrsmStrPtrVec2D decoy4_proteoforms;
  {
    PrsmStrPtrVec* targets[] = {&target_prsms1, &target_prsms2, &target_prsms3,
                                &target_prsms4};
    PrsmStrPtrVec* decoys[] = {&decoy_prsms1, &decoy_prsms2, &decoy_prsms3,
                               &decoy_prsms4};
    PrsmStrPtrVec* alls[] = {&all_prsms1, &all_prsms2, &all_prsms3,
                             &all_prsms4};
    PrsmStrPtrVec2D* target_forms[] = {&target1_proteoforms,
                                       &target2_proteoforms,
                                       &target3_proteoforms,
                                       &target4_proteoforms};
    PrsmStrPtrVec2D* decoy_forms[] = {&decoy1_proteoforms, &decoy2_proteoforms,
                                      &decoy3_proteoforms, &decoy4_proteoforms};
    for (int g = 0; g < 4; g++) {
      *alls[g] = *targets[g];
      alls[g]->insert(alls[g]->end(), decoys[g]->begin(), decoys[g]->end());
      std::sort(alls[g]->begin(), alls[g]->end(),
                PrsmStr::cmpEValueIncProtInc);
      PrsmStrPtrVec2D all_proteoforms = prsm_fdr::getProteoformGroups(*alls[g]);
      prsm_fdr::splitTargetDecoyProteoforms(all_proteoforms, *target_forms[g],
                                            *decoy_forms[g]);
    }
  }

  prsm_fdr::computeProteoformFdr(target1_proteoforms, decoy1_proteoforms);
  prsm_fdr::computeProteoformFdr(target2_proteoforms, decoy2_proteoforms);
  prsm_fdr::computeProteoformFdr(target3_proteoforms, decoy3_proteoforms);
  prsm_fdr::computeProteoformFdr(target4_proteoforms, decoy4_proteoforms);

  // protein-level FDR from the combined proteoforms of all four groups
  PrsmStrPtrVec2D target_proteoforms_all;
  PrsmStrPtrVec2D decoy_proteoforms_all;
  for (auto* v : {&target1_proteoforms, &target2_proteoforms,
                  &target3_proteoforms, &target4_proteoforms}) {
    target_proteoforms_all.insert(target_proteoforms_all.end(), v->begin(),
                                  v->end());
  }
  for (auto* v : {&decoy1_proteoforms, &decoy2_proteoforms,
                  &decoy3_proteoforms, &decoy4_proteoforms}) {
    decoy_proteoforms_all.insert(decoy_proteoforms_all.end(), v->begin(),
                                 v->end());
  }
  auto cmp_form = [](const PrsmStrPtrVec& a, const PrsmStrPtrVec& b) {
    return PrsmStr::cmpEValueIncProtInc(a[0], b[0]);
  };
  std::sort(target_proteoforms_all.begin(), target_proteoforms_all.end(),
            cmp_form);
  std::sort(decoy_proteoforms_all.begin(), decoy_proteoforms_all.end(),
            cmp_form);
  PrsmStrPtrVec2D target_proteins =
      prsm_fdr::getProteinGroups(target_proteoforms_all);
  PrsmStrPtrVec2D decoy_proteins =
      prsm_fdr::getProteinGroups(decoy_proteoforms_all);
  prsm_fdr::computeProteinFdr(target_proteins, decoy_proteins);

  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);

  // every PrSM, targets and decoys, group by group in spectrum order
  PrsmStrPtrVec* alls[] = {&all_prsms1, &all_prsms2, &all_prsms3,
                           &all_prsms4};
  for (int g = 0; g < 4; g++) {
    std::sort(alls[g]->begin(), alls[g]->end(),
              PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
    writer.writeVector(*alls[g]);
  }
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
