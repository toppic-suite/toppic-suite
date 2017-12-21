//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include <string>
#include <algorithm>
#include <vector>

#include "base/file_util.hpp"
#include "base/proteoform_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_species.hpp"

namespace prot {

std::vector<PrsmStrPtrVec> PrsmSpecies::groupProteins(const PrsmStrPtrVec &prsm_ptrs) {
  // get max shift number
  int max_shift_number = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int cur_shift_number = prsm_ptrs[i]->getUnexpectedPtmNum();
    if (max_shift_number < cur_shift_number) {
      max_shift_number = cur_shift_number;
    }
  }
  // get proteoform groups
  std::vector<PrsmStrPtrVec> proteogroups;
  for (int shift = 0; shift <= max_shift_number; shift++) {
    PrsmStrPtrVec proteo_ptrs;
    for (size_t i = 0; i < prsm_ptrs.size(); i++) {
      if (shift == prsm_ptrs[i]->getUnexpectedPtmNum()) {
        proteo_ptrs.push_back(prsm_ptrs[i]);
      }
    }
    proteogroups.push_back(proteo_ptrs);
  }
  return proteogroups;
}

std::vector<PrsmStrPtrVec> PrsmSpecies::getZeroPtmList(const PrsmStrPtrVec& proteo_ptrs, double ppo) {
  std::vector<PrsmStrPtrVec> species;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    bool is_found = false;
    for (size_t j = 0; j < species.size(); j++) {
      if (PrsmStr::isSameSeqAndMass(proteo_ptrs[i], species[j][0], ppo)) {
        species[j].push_back(proteo_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmStrPtrVec new_species;
      new_species.push_back(proteo_ptrs[i]);
      species.push_back(new_species);
    }
  }
  return species;
}

void PrsmSpecies::setProtId(PrsmStrPtrVec& prsm_ptrs) {
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

void PrsmSpecies::setSpeciesId(PrsmStrPtrVec& prsm_ptrs, double ppo) {
  std::vector<PrsmStrPtrVec> proteo_groups = groupProteins(prsm_ptrs);

  // find zero ptm species
  std::vector<PrsmStrPtrVec> species = getZeroPtmList(proteo_groups[0], ppo);

  for (size_t i = 1; i < proteo_groups.size(); i++) {
    for (size_t j = 0; j < proteo_groups[i].size(); j++) {
      bool is_found = false;
      for (size_t m = 0; m < species.size(); m++) {
        if (PrsmStr::isStrictCompatiablePtmSpecies(proteo_groups[i][j], species[m][0], ppo)) {
          species[m].push_back(proteo_groups[i][j]);
          is_found = true;
          break;
        }
      }
      if (!is_found) {
        PrsmStrPtrVec new_species;
        new_species.push_back(proteo_groups[i][j]);
        species.push_back(new_species);
      }
    }
  }

  for (size_t i = 0; i < species.size(); i++) {
    for (size_t j = 0; j < species[i].size(); j++) {
      species[i][j]->setClusterId(i);
    }
  }
}

void PrsmSpecies::process() {
  std::string base_name = file_util::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmStrPtrVec prsm_ptrs = PrsmReader::readAllPrsmStrs(input_file_name);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpEValueInc);
  setProtId(prsm_ptrs);
  setSpeciesId(prsm_ptrs, ppo_);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpSpectrumIdIncPrecursorIdInc);
  // output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prot


