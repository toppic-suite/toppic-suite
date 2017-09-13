// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <string>
#include <algorithm>
#include <vector>

#include "base/file_util.hpp"
#include "base/proteoform_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_species.hpp"

namespace prot {

ProteoformPtrVec2D PrsmSpecies::groupProteins(const PrsmPtrVec &prsm_ptrs) {
  // get max shift number
  int max_shift_number = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int cur_shift_number
        = prsm_ptrs[i]->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED);
    if (max_shift_number < cur_shift_number) {
      max_shift_number = cur_shift_number;
    }
  }
  // get proteoform groups
  ProteoformPtrVec2D proteogroups;
  for (int shift = 0; shift <= max_shift_number; shift++) {
    ProteoformPtrVec proteo_ptrs;
    for (size_t i = 0; i < prsm_ptrs.size(); i++) {
      if (shift == prsm_ptrs[i]->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED)) {
        proteo_ptrs.push_back(prsm_ptrs[i]->getProteoformPtr());
      }
    }
    proteogroups.push_back(proteo_ptrs);
  }
  return proteogroups;
}

ProteoformPtrVec2D PrsmSpecies::getZeroPtmList(const ProteoformPtrVec& proteo_ptrs, double ppo) {
  ProteoformPtrVec2D species;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    bool is_found = false;
    for (size_t j = 0; j < species.size(); j++) {
      if (ProteoformUtil::isSameSeqAndMass(proteo_ptrs[i], species[j][0], ppo)) {
        species[j].push_back(proteo_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      ProteoformPtrVec new_species;
      new_species.push_back(proteo_ptrs[i]);
      species.push_back(new_species);
    }
  }
  return species;
}

void PrsmSpecies::setProtId(PrsmPtrVec& prsm_ptrs) {
  PrsmPtrVec2D proteins;
  std::vector<std::string> protein_names;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    std::string name = prsm_ptrs[i]->getProteoformPtr()->getSeqName();
    bool is_found = false;
    for (size_t j = 0; j < protein_names.size(); j++) {
      if (protein_names[j] == name) {
        proteins[j].push_back(prsm_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmPtrVec new_protein;
      new_protein.push_back(prsm_ptrs[i]);
      proteins.push_back(new_protein);
      protein_names.push_back(name);
    }
  }

  for (size_t i = 0; i < proteins.size(); i++) {
    for (size_t j = 0; j < proteins[i].size(); j++) {
      proteins[i][j]->getProteoformPtr()->setProtId(i);
    }
  }
}

void PrsmSpecies::setSpeciesId(PrsmPtrVec& prsm_ptrs, double ppo) {
  ProteoformPtrVec2D proteo_groups = groupProteins(prsm_ptrs);

  // find zero ptm species
  ProteoformPtrVec2D species = getZeroPtmList(proteo_groups[0], ppo);

  for (size_t i = 1; i < proteo_groups.size(); i++) {
    for (size_t j = 0; j < proteo_groups[i].size(); j++) {
      bool is_found = false;
      for (size_t m = 0; m< species.size(); m++) {
        if (ProteoformUtil::isStrictCompatiablePtmSpecies(proteo_groups[i][j], species[m][0], ppo)) {
          species[m].push_back(proteo_groups[i][j]);
          is_found = true;
          break;
        }
      }
      if (!is_found) {
        ProteoformPtrVec new_species;
        new_species.push_back(proteo_groups[i][j]);
        species.push_back(new_species);
      }
    }
  }

  for (size_t i = 0; i < species.size(); i++) {
    for (size_t j = 0; j < species[i].size(); j++) {
      species[i][j]->setSpeciesId(i);
    }
  }
}

void PrsmSpecies::process() {
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;

  PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name_,
                                                  fix_mod_ptr_vec_);
  // sort(prsm_ptrs.begin(),prsm_ptrs.end(),Prsm::cmpSpectrumIdIncPrecursorIdInc);
  // sort(prsm_ptrs.begin(),prsm_ptrs.end(),Prsm::cmpMatchFragmentDecMatchPeakDec);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), Prsm::cmpEValueInc);
  setProtId(prsm_ptrs);
  setSpeciesId(prsm_ptrs, ppo_);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), Prsm::cmpSpectrumIdIncPrecursorIdInc);
  // output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}


}  // namespace prot


