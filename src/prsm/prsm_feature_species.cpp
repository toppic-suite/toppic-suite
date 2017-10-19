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

#include <boost/algorithm/string.hpp>

#include "base/file_util.hpp"
#include "base/proteoform_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_feature_species.hpp"
#include "spec/msalign_reader.hpp"

namespace prot {

void PrsmFeatureSpecies::setProtId(PrsmStrPtrVec& prsm_ptrs) {
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

void PrsmFeatureSpecies::setSpeciesId(PrsmStrPtrVec& prsm_ptrs) {
  std::vector<PrsmStrPtrVec> species;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    bool is_found = false;
    PrsmStrPtr cur_ptr = prsm_ptrs[i];
    for (size_t j = 0; j < species.size(); j++) {
      PrsmStrPtr ref_ptr = species[j][0];
      if (cur_ptr->getProtId() == ref_ptr->getProtId()) {
        if (cur_ptr->getPrecFeatureId() == ref_ptr->getPrecFeatureId()) {
          species[j].push_back(cur_ptr);
          is_found = true;
          break;
        }
        if (std::abs(cur_ptr->getAdjustedPrecMass() - ref_ptr->getAdjustedPrecMass()) <= prec_error_tole_) {
          species[j].push_back(cur_ptr);
          is_found = true;
          break;
        }
      } else if (cur_ptr->getProteinMatchSeq() == ref_ptr->getProteinMatchSeq()) {
        species[j].push_back(cur_ptr);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      PrsmStrPtrVec new_species;
      new_species.push_back(prsm_ptrs[i]);
      species.push_back(new_species);
    }
  }
  for (size_t i = 0; i < species.size(); i++) {
    for (size_t j = 0; j < species[i].size(); j++) {
      species[i][j]->setSpeciesId(i);
    }
  }
}

void PrsmFeatureSpecies::process() {
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmStrPtrVec prsm_ptrs = PrsmReader::readAllPrsmStrs(input_file_name);
  PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(input_file_name);
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name_);

  // read TopFD featuers
  std::vector<int> feature_spec_ids;
  std::vector<int> feature_ids;
  std::vector<double> feature_intens;
  std::ifstream infile(feature_file_name_);
  std::string line;
  while (std::getline(infile, line)) {
    if (line[0] == '#' || line == "" || line[0] == 'I') {
      continue;
    }
    std::vector<std::string> strs; 
    boost::split(strs, line, boost::is_any_of("\t "));
    feature_spec_ids.push_back(std::stoi(strs[0]));
    feature_ids.push_back(std::stoi(strs[6]));
    feature_intens.push_back(std::stod(strs[7]));
  }
  infile.close();

  size_t k = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int spec_id = prsm_ptrs[i]->getSpectrumId();
    while (feature_spec_ids[k] != spec_id) {k++;}
    prsm_ptrs[i]->setPrecFeatureId(feature_ids[k]);
    prsm_ptrs[i]->setPrecFeatureInte(feature_intens[k]);

    PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(seq_reader, fix_mod_ptr_vec_);
    if (prsm_ptr != nullptr) {
      prsm_ptrs[i]->setProteinMatchSeq(prsm_ptr->getProteoformPtr()->getProteinMatchSeq());
    }
  }
  prsm_reader->close();
  prsm_reader = nullptr;
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpEValueInc);
  setProtId(prsm_ptrs);
  setSpeciesId(prsm_ptrs);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpSpectrumIdIncPrecursorIdInc);
  // output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prot
