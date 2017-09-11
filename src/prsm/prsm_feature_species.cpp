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


#include "base/file_util.hpp"
#include "base/proteoform_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_feature_species.hpp"
#include "spec/msalign_reader.hpp"

namespace prot {

void PrsmFeatureSpecies::setProtId(PrsmPtrVec& prsm_ptrs){
  PrsmPtrVec2D proteins;
  std::vector<std::string> protein_names;
  for(size_t i = 0; i < prsm_ptrs.size(); i++) {
    std::string name = prsm_ptrs[i]->getProteoformPtr()->getSeqName();
    bool is_found = false;
    for(size_t j = 0; j < protein_names.size(); j++){
      if(protein_names[j] == name) {
        proteins[j].push_back(prsm_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if(!is_found){
      PrsmPtrVec new_protein;
      new_protein.push_back(prsm_ptrs[i]);
      proteins.push_back(new_protein);
      protein_names.push_back(name);
    }
  }

  for(size_t i = 0; i < proteins.size(); i++){
    for (size_t j = 0; j < proteins[i].size(); j++) {
      proteins[i][j]->getProteoformPtr()->setProtId(i);
    }
  }
}

void PrsmFeatureSpecies::setSpeciesId(const PrsmPtrVec& prsm_ptrs){
  PrsmPtrVec2D species;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    bool is_found = false;
    PrsmPtr cur_ptr = prsm_ptrs[i];
    for (size_t j = 0; j < species.size(); j++) {
      PrsmPtr ref_ptr = species[j][0];
      if (cur_ptr->getProteoformPtr()->getProtId() == ref_ptr->getProteoformPtr()->getProtId()) {
        if (cur_ptr->getPrecFeatureId() == ref_ptr->getPrecFeatureId()) {
          species[j].push_back(cur_ptr);
          is_found = true;
          break;
        }
        if (std::abs (cur_ptr->getAdjustedPrecMass() - ref_ptr->getAdjustedPrecMass()) <= prec_error_tole_) {
          species[j].push_back(cur_ptr);
          is_found = true;
          break;
        }
      }
    }
    if (!is_found) {
      PrsmPtrVec new_species;
      new_species.push_back(prsm_ptrs[i]);
      species.push_back(new_species);
    }
  }
  for(size_t i = 0; i < species.size(); i++){
    for (size_t j = 0; j < species[i].size(); j++) {
      species[i][j]->getProteoformPtr()->setSpeciesId(i);
    }
  }
}

void PrsmFeatureSpecies::process() {
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name_,
                                                  fix_mod_ptr_vec_);
  MsAlignReaderPtr msreader
      = std::make_shared<MsAlignReader>(feature_file_name_,
                                        prsm_para_ptr_->getGroupSpecNum(),
                                        prsm_para_ptr_->getSpParaPtr()->getActivationPtr());
  SpectrumSetPtr spec_set_ptr = msreader->getNextSpectrumSet(prsm_para_ptr_->getSpParaPtr());
  LOG_DEBUG("prsm_ptrs.size() " << prsm_ptrs.size());
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int spec_id = prsm_ptrs[i]->getSpectrumId();
    while (spec_set_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr()->getId() != spec_id) {
      spec_set_ptr = msreader->getNextSpectrumSet(prsm_para_ptr_->getSpParaPtr()); 
    }
    prsm_ptrs[i]->setPrecFeatureId(spec_set_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr()->getFeatureId());
    prsm_ptrs[i]->setPrecFeatureInte(spec_set_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr()->getFeatureInte());
  } 
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), Prsm::cmpEValueInc);
  setProtId(prsm_ptrs);
  setSpeciesId(prsm_ptrs);
  sort(prsm_ptrs.begin(), prsm_ptrs.end(), Prsm::cmpSpectrumIdIncPrecursorIdInc);
  //output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prot
