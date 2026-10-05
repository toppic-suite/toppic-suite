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

#include "search/graph/proteo_anno.hpp"

#include <memory>
#include <string>
#include <vector>

#include "common/base/amino_acid_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/prot_mod_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/residue_util.hpp"
#include "common/util/logger.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "seq/alter_type.hpp"

namespace toppic {

ProteoAnno::ProteoAnno(const ModPtrVec& fix_mod_ptr_vec,
                       const ProtModPtrVec& prot_mod_ptr_vec,
                       const ModPtrVec& var_mod_ptr_vec) {
  fix_mod_ptr_vec_ = fix_mod_ptr_vec;
  prot_mod_ptr_vec_ = prot_mod_ptr_vec;
  var_mod_ptr_vec_ = var_mod_ptr_vec;
  is_nme_ = false;
  for (size_t i = 0; i < var_mod_ptr_vec_.size(); i++) {
    ResiduePtr res_ptr = var_mod_ptr_vec[i]->getModResiduePtr();
    AminoAcidPtr acid_ptr = res_ptr->getAminoAcidPtr();
    if (ptm_map_.find(acid_ptr) == ptm_map_.end()) {
      ResiduePtrVec cur_vec;
      cur_vec.push_back(res_ptr);
      ptm_map_[acid_ptr] = cur_vec;
    } else {
      ptm_map_[acid_ptr].push_back(res_ptr);
    }
  }
}

// The residues come from the sequence's residue/PTM pairs: a PTM written in
// the FASTA sequence (M[Oxidation]) makes the modified residue, as
// proteoform_factory::geneDbProteoformPtr makes it for the database
// proteoform. The raw string used to be read instead, whose brackets were
// dropped and whose PTM name became residues (M[Oxidation] as M and O).
void ProteoAnno::anno(const StringPairVec& acid_ptm_pairs, bool is_complete) {
  ResiduePtrVec residue_ptr_vec =
      residue_util::convertStrToResiduePtrVec(acid_ptm_pairs, fix_mod_ptr_vec_);
  res_vec_2d_.clear();
  shift_vec_2d_.clear();
  // A residue with a PTM, of a fixed modification or of the FASTA sequence,
  // is FIXED: its one edge, which getMassGraphPtr keeps alone, so that it
  // takes no variable PTM or protein N-terminal form. An unmodified residue
  // is marked INPUT.
  for (size_t i = 0; i < residue_ptr_vec.size(); i++) {
    ResiduePtrVec cur_res_vec;
    std::vector<int> cur_shift_vec;
    ResiduePtr res_ptr = residue_ptr_vec[i];
    cur_res_vec.push_back(res_ptr);
    if (PtmBase::isEmptyPtmPtr(res_ptr->getPtmPtr())) {
      cur_shift_vec.push_back(AlterType::INPUT->getId());
    } else {
      cur_shift_vec.push_back(AlterType::FIXED->getId());
    }
    res_vec_2d_.push_back(cur_res_vec);
    shift_vec_2d_.push_back(cur_shift_vec);
  }
  LOG_DEBUG("input complete");

  // protein N-terminal forms, of the protein's N-terminal sequence only and
  // those allowMod accepts for it. Their edges are at residue 0, so that
  // only an alignment that starts at the protein N-terminus (vertex 0) can
  // use them; GraphPostProcessor reads the form back from the first
  // variable shift and the start. NME is the start at vertex 1, and is_nme_
  // tells whether that start is the protein N-terminus.
  // A first residue with a PTM has no form edge (its edge is FIXED): the
  // acetylated forms are not open to it, NME (a start, not an edge) is.
  is_nme_ = false;
  bool first_fixed = !residue_ptr_vec.empty() &&
                     !PtmBase::isEmptyPtmPtr(residue_ptr_vec[0]->getPtmPtr());
  for (size_t i = 0; is_complete && i < prot_mod_ptr_vec_.size(); i++) {
    ProtModPtr mod_ptr = prot_mod_ptr_vec_[i];
    if (!prot_mod_util::allowMod(mod_ptr, residue_ptr_vec)) {
      continue;
    }
    LOG_DEBUG("i " << i << " mod " << mod_ptr);
    if (mod_ptr->getType() == ProtModBase::getType_NME()) {
      is_nme_ = true;
    } else if (first_fixed) {
      continue;
    } else if (mod_ptr->getType() == ProtModBase::getType_M_ACETYLATION()) {
      // the acetylated M
      ResiduePtr mut_residue_ptr = mod_ptr->getModPtr()->getModResiduePtr();
      res_vec_2d_[0].push_back(mut_residue_ptr);
      shift_vec_2d_[0].push_back(AlterType::PROTEIN_VARIABLE->getId());
    } else if (mod_ptr->getType() == ProtModBase::getType_NME_ACETYLATION()) {
      // the M removed and residue 1 acetylated: an edge in place of the M
      // that weighs the acetyl group alone. It used to be an acetylated
      // residue 1, which an alignment keeping the M could use too.
      PtmPtr ptm_ptr = mod_ptr->getModPtr()->getModResiduePtr()->getPtmPtr();
      ResiduePtr acetyl_ptr = std::make_shared<Residue>(
          AminoAcidBase::getEmptyAminoAcidPtr(), ptm_ptr);
      res_vec_2d_[0].push_back(acetyl_ptr);
      shift_vec_2d_[0].push_back(AlterType::PROTEIN_VARIABLE->getId());
    }
    LOG_DEBUG("round complete");
  }
  LOG_DEBUG("protein mod complete");

  // variable ptms, of the residues without a PTM
  for (size_t i = 0; i < residue_ptr_vec.size(); i++) {
    if (!PtmBase::isEmptyPtmPtr(residue_ptr_vec[i]->getPtmPtr())) {
      continue;
    }
    AminoAcidPtr acid_ptr = residue_ptr_vec[i]->getAminoAcidPtr();
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
          shift_vec_2d_[i].push_back(AlterType::VARIABLE->getId());
        }
      }
    }
  }
  LOG_DEBUG("variable mod complete");
}
}  // namespace toppic
