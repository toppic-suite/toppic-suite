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
#include <vector>

#include "base/xml_dom.hpp"
#include "base/change_type.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"
#include "base/ptm_base.hpp"
#include "base/residue_base.hpp"
#include "base/residue_util.hpp"
#include "base/xml_dom_document.hpp"
#include "graph/proteo_anno.hpp"

namespace prot {

ProteoAnno::ProteoAnno(const ModPtrVec &fix_mod_ptr_vec,
                       const ProtModPtrVec &prot_mod_ptr_vec,
                       const ModPtrVec &var_mod_ptr_vec) {
  fix_mod_ptr_vec_ = fix_mod_ptr_vec;
  prot_mod_ptr_vec_ = prot_mod_ptr_vec;
  var_mod_ptr_vec_ = var_mod_ptr_vec;
  is_nme_ = false;
  for (size_t i = 0; i < var_mod_ptr_vec_.size(); i++) {
    ResiduePtr res_ptr = var_mod_ptr_vec[i]->getModResiduePtr();
    AcidPtr acid_ptr = res_ptr->getAcidPtr();
    if (ptm_map_.find(acid_ptr) == ptm_map_.end()) {
      ResiduePtrVec cur_vec;
      cur_vec.push_back(res_ptr);
      ptm_map_[acid_ptr]= cur_vec;
    } else {
      ptm_map_[acid_ptr].push_back(res_ptr);
    }
  }
}

void ProteoAnno::anno(const std::string &seq, bool is_complete) {
  ResiduePtrVec residue_ptr_vec = ResidueUtil::convertStrToResiduePtrVec(seq, fix_mod_ptr_vec_);
  res_vec_2d_.clear();
  change_vec_2d_.clear();
  // input and fixed mod
  for (size_t i = 0; i < residue_ptr_vec.size(); i++) {
    ResiduePtrVec cur_res_vec;
    std::vector<int> cur_change_vec;
    ResiduePtr res_ptr = residue_ptr_vec[i];
    cur_res_vec.push_back(res_ptr);
    if (PtmBase::isEmptyPtmPtr(res_ptr->getPtmPtr())) {
      cur_change_vec.push_back(ChangeType::INPUT->getId());
    } else {
      cur_change_vec.push_back(ChangeType::FIXED->getId());
    }
    res_vec_2d_.push_back(cur_res_vec);
    change_vec_2d_.push_back(cur_change_vec);
  }
  LOG_DEBUG("input complete");

  // protein mod
  for (size_t i = 0; i < prot_mod_ptr_vec_.size(); i++) {
    ProtModPtr mod_ptr = prot_mod_ptr_vec_[i];
    if (!ProtModUtil::allowMod(mod_ptr, residue_ptr_vec)) {
      continue;
    }
    LOG_DEBUG("i " << i << " mod " << mod_ptr);
    if (is_complete && mod_ptr->getType() == ProtModBase::getType_NME()) {
      LOG_DEBUG("NME");
      // add empty residue to the first methinine residue
      is_nme_ = true;
      ResiduePtr empty_residue_ptr = ResidueBase::getEmptyResiduePtr();
      // LOG_DEBUG("empty acid mass " << acid_ptr->getMonoMass());
      // LOG_DEBUG("empty ptm mass " << ptm_ptr->getMonoMass());
      LOG_DEBUG("empty residue mass " << empty_residue_ptr->getMass());
      if (empty_residue_ptr == nullptr) {
        LOG_ERROR("Proteoform:: residue not found");
        throw("Residue not found");
      }
      res_vec_2d_[0].push_back(empty_residue_ptr);
      change_vec_2d_[0].push_back(ChangeType::PROTEIN_VARIABLE->getId());
    } else if (is_complete && mod_ptr->getType() == ProtModBase::getType_M_ACETYLATION()) {
      ResiduePtr mut_residue_ptr = mod_ptr->getModPtr()->getModResiduePtr();
      res_vec_2d_[0].push_back(mut_residue_ptr);
      change_vec_2d_[0].push_back(ChangeType::PROTEIN_VARIABLE->getId());
    } else if (is_complete && mod_ptr->getType() == ProtModBase::getType_NME_ACETYLATION()) {
      LOG_DEBUG("NME_ACETYLATION");
      // add acetylation to the second residue
      is_nme_ = true;
      ResiduePtr mut_residue_ptr = mod_ptr->getModPtr()->getModResiduePtr();
      res_vec_2d_[1].push_back(mut_residue_ptr);
      change_vec_2d_[1].push_back(ChangeType::PROTEIN_VARIABLE->getId());
    }
    LOG_DEBUG("round complete");
  }
  LOG_DEBUG("protein mod complete");

  // variable ptms
  for (size_t i = 0; i < residue_ptr_vec.size(); i++) {
    AcidPtr acid_ptr = residue_ptr_vec[i]->getAcidPtr();
    // if exist modified residues
    if (ptm_map_.count(acid_ptr)) {
      ResiduePtrVec mod_res_vec = ptm_map_[acid_ptr];
      // remove duplications
      ResiduePtrVec exist_res_vec = res_vec_2d_[i];
      for (size_t j = 0; j < mod_res_vec.size(); j++) {
        bool found = false;
        for (size_t k = 0; k < exist_res_vec.size(); k++) {
          if (mod_res_vec[j] == exist_res_vec[k]) {
            found = true;
            break;
          }
        }
        if (!found) {
          res_vec_2d_[i].push_back(mod_res_vec[j]);
          change_vec_2d_[i].push_back(ChangeType::VARIABLE->getId());
        }
      }
    }
  }
  LOG_DEBUG("variable mod complete");
}
}  // namespace prot
