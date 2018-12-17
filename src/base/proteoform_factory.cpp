//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#include <sstream>
#include <algorithm>
#include <string>

#include "base/logger.hpp"
#include "base/fasta_reader.hpp"
#include "base/ptm_base.hpp"
#include "base/residue_base.hpp"
#include "base/residue_util.hpp"
#include "base/residue_seq.hpp"
#include "base/trunc_util.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"
#include "base/proteoform_factory.hpp"

namespace toppic {

namespace proteoform_factory {

ProteoformPtr geneDbProteoformPtr(FastaSeqPtr fasta_seq_ptr, ModPtrVec fix_mod_list) {
  if (fasta_seq_ptr == nullptr) {
    return ProteoformPtr(nullptr);
  }
  ProtModPtr none_prot_mod_ptr = ProtModBase::getProtModPtr_NONE();
  ResiduePtrVec residue_ptrs = residue_util::convertStrToResiduePtrVec(fasta_seq_ptr->getAcidPtmPairVec());
  int start_pos = 0;
  int end_pos = static_cast<int>(residue_ptrs.size()) - 1;

  MassShiftPtrVec shift_list;
  // add input ptms;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    if (residue_ptrs[i]->getPtmPtr() != PtmBase::getEmptyPtmPtr()) {
      ResiduePtr ori_residue = ResidueBase::getBaseResiduePtr(residue_ptrs[i]->getAminoAcidPtr());
      ModPtr mod_ptr = ModBase::getBaseModPtr(ori_residue, residue_ptrs[i]);
      ChangePtr change_ptr
          = std::make_shared<Change>(i, i + 1, MassShiftType::INPUT, mod_ptr->getShift(), mod_ptr);
      MassShiftPtr shift_ptr
          = std::make_shared<MassShift>(i, i + 1, change_ptr->getTypePtr());
      shift_ptr->setChangePtr(change_ptr);
      shift_list.push_back(shift_ptr);
    }
  }

  // add fixed ptms;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    for (size_t j = 0; j < fix_mod_list.size(); j++) {
      if (residue_ptrs[i] == fix_mod_list[j]->getOriResiduePtr()) {
        residue_ptrs[i] = fix_mod_list[j]->getModResiduePtr();
        ChangePtr change_ptr
            = std::make_shared<Change>(i, i + 1, MassShiftType::FIXED, fix_mod_list[j]->getShift(), fix_mod_list[j]);
        MassShiftPtr shift_ptr
            = std::make_shared<MassShift>(i, i + 1, change_ptr->getTypePtr());
        shift_ptr->setChangePtr(change_ptr);
        shift_list.push_back(shift_ptr);
        break;
      }
    }
  }

  ResSeqPtr res_seq_ptr = std::make_shared<ResidueSeq>(residue_ptrs);
  return std::make_shared<Proteoform>(fasta_seq_ptr, none_prot_mod_ptr, start_pos,
                                      end_pos, res_seq_ptr, shift_list);
}

ProteoformPtr geneProtModProteoform(ProteoformPtr db_form_ptr, ProtModPtr prot_mod_ptr) {
  // check if the proteoform can be truncated
  ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
  bool valid_mod = prot_mod_util::allowMod(prot_mod_ptr, db_res_seq_ptr->getResidues());
  if (!valid_mod) {
    // LOG_DEBUG("NO valid mod");
    return ProteoformPtr(nullptr);
  }

  TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
  int start = trunc_ptr->getTruncLen();
  /* last bp index */
  int end = db_form_ptr->getLen();
  // copy input changes
  MassShiftPtrVec ori_shift_ptrs = db_form_ptr->getMassShiftPtrVec();
  MassShiftPtrVec shift_ptrs;
  for (size_t i = 0; i < ori_shift_ptrs.size(); i++) {
    if (ori_shift_ptrs[i]->getLeftBpPos() >= start
        && ori_shift_ptrs[i]->getRightBpPos() <= end + 1) {
      MassShiftPtr shift_ptr = std::make_shared<MassShift>(ori_shift_ptrs[i]->getLeftBpPos(),
                                                           ori_shift_ptrs[i]->getRightBpPos(),
                                                           ori_shift_ptrs[i]->getTypePtr());
      ChangePtrVec change_ptrs = ori_shift_ptrs[i]->getChangePtrVec();
      for (size_t k = 0; k < change_ptrs.size(); k++) {
        ChangePtr change_ptr = Change::geneChangePtr(change_ptrs[k], start); 
        shift_ptr->setChangePtr(change_ptr);
      }
      shift_ptrs.push_back(shift_ptr);
    } 
  }

  ResiduePtrVec ori_vec = db_res_seq_ptr->getResidues();
  ResiduePtrVec new_vec(ori_vec.begin()+start, ori_vec.begin()+end);
  // apply mod
  ModPtr mod_ptr = prot_mod_ptr->getModPtr();
  if (!ModBase::isNoneModPtr(mod_ptr)) {
    int mod_pos = prot_mod_ptr->getModPos();
    int new_pos = mod_pos - start;
    new_vec[new_pos] = prot_mod_ptr->getModPtr()->getModResiduePtr();
    ChangePtr change_ptr
        = std::make_shared<Change>(new_pos, new_pos + 1, MassShiftType::PROTEIN_VARIABLE,
                                   mod_ptr->getShift(), mod_ptr);
    MassShiftPtr shift_ptr = std::make_shared<MassShift>(change_ptr->getLeftBpPos(),
                                                         change_ptr->getRightBpPos(),
                                                         change_ptr->getTypePtr());
    shift_ptr->setChangePtr(change_ptr);
    shift_ptrs.push_back(shift_ptr);
  }
  ResSeqPtr seq_ptr = std::make_shared<ResidueSeq>(new_vec);

  FastaSeqPtr fasta_seq_ptr = db_form_ptr->getFastaSeqPtr();

  return std::make_shared<Proteoform>(db_form_ptr->getFastaSeqPtr(), prot_mod_ptr, start,
                                      db_res_seq_ptr->getLen() - 1, seq_ptr, shift_ptrs);
}

ProteoformPtr geneSubProteoform(ProteoformPtr proteoform_ptr, int local_start, int local_end) {
  ResiduePtrVec residues;
  ResSeqPtr res_seq_ptr = proteoform_ptr->getResSeqPtr();
  for (int i = local_start; i <= local_end; i++) {
    residues.push_back(res_seq_ptr->getResiduePtr(i));
  }
  ResSeqPtr seq_ptr = std::make_shared<ResidueSeq>(residues);
  MassShiftPtrVec shift_list;
  MassShiftPtrVec ori_shift_list = proteoform_ptr->getMassShiftPtrVec();

  for (size_t i = 0; i < ori_shift_list.size(); i++) {
    if (ori_shift_list[i]->getLeftBpPos() >= local_start
        && ori_shift_list[i]->getRightBpPos() <= local_end + 1) {
      MassShiftPtr shift_ptr = std::make_shared<MassShift>(ori_shift_list[i]->getLeftBpPos(),
                                                           ori_shift_list[i]->getRightBpPos(),
                                                           ori_shift_list[i]->getTypePtr());
      ChangePtrVec change_ptrs = ori_shift_list[i]->getChangePtrVec();
      for (size_t k = 0; k < change_ptrs.size(); k++) {
        ChangePtr change_ptr = Change::geneChangePtr(change_ptrs[k], local_start); 
        shift_ptr->setChangePtr(change_ptr);
      }
      shift_list.push_back(shift_ptr);
    } 
  }

  ProtModPtr prot_mod_ptr = proteoform_ptr->getProtModPtr();

  return std::make_shared<Proteoform>(proteoform_ptr->getFastaSeqPtr(), prot_mod_ptr,
                                      local_start + proteoform_ptr->getStartPos(),
                                      local_end + proteoform_ptr->getStartPos(),
                                      seq_ptr, shift_list);
}

ProteoformPtr geneProteoform(ProteoformPtr proteoform, int start_pos, int end_pos,
                             const MassShiftPtrVec & mass_shift_vec,
                             const ModPtrVec & mod_ptr_vec) {
  ResSeqPtr new_res_seq = std::make_shared<ResidueSeq>(
      residue_util::convertStrToResiduePtrVec(
          proteoform->getFastaSeqPtr()->getRawSeq().substr(start_pos, end_pos - start_pos + 1), 
          mod_ptr_vec));
  for (size_t i = 0; i < mass_shift_vec.size(); i++) {
    int left = proteoform->getStartPos() + mass_shift_vec[i]->getLeftBpPos() - start_pos;
    int right = proteoform->getStartPos() + mass_shift_vec[i]->getRightBpPos() - start_pos;
    mass_shift_vec[i]->setLeftBpPos(left);
    mass_shift_vec[i]->setRightBpPos(right);
  }

  return std::make_shared<Proteoform>(proteoform->getFastaSeqPtr(),
                                      proteoform->getProtModPtr(), start_pos, end_pos, 
                                      new_res_seq, mass_shift_vec);
}

ProteoformPtrVec geneProtModProteoform(ProteoformPtr proteo_ptr, const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (size_t j = 0; j < prot_mods.size(); j++) {
    ProteoformPtr ptr = geneProtModProteoform(proteo_ptr, prot_mods[j]);
    if (ptr.get() != nullptr) {
      new_forms.push_back(ptr);
    }
  }
  return new_forms;
}

ProteoformPtrVec geneProtModProteoform(const ProteoformPtrVec &ori_forms,
                                       const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (size_t i = 0; i < ori_forms.size(); i++) {
    for (size_t j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = geneProtModProteoform(ori_forms[i], prot_mods[j]);
      if (ptr.get() != nullptr) {
        new_forms.push_back(ptr);
      }
    }
  }
  return new_forms;
}

ProteoformPtrVec2D gene2DProtModProteoform(const ProteoformPtrVec &ori_forms,
                                           const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec2D new_forms;
  for (size_t i = 0; i < ori_forms.size(); i++) {
    ProteoformPtrVec mod_forms;
    for (size_t j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = geneProtModProteoform(ori_forms[i], prot_mods[j]);
      if (ptr.get() != nullptr) {
        mod_forms.push_back(ptr);
      }
    }
    new_forms.push_back(mod_forms);
  }
  return new_forms;
}

ProteoformPtrVec readFastaToProteoformPtrVec(const std::string &file_name,
                                             const ModPtrVec &fix_mod_list) {
  LOG_DEBUG("start open file " << file_name);
  FastaReader reader(file_name);
  LOG_DEBUG("open file done " << file_name);

  ProteoformPtrVec list;
  FastaSeqPtr seq_ptr = reader.getNextSeq();
  int count = 0;
  while (seq_ptr != nullptr) {
    ProteoformPtr proteo_ptr = geneDbProteoformPtr(seq_ptr, fix_mod_list);
    list.push_back(proteo_ptr);
    seq_ptr = reader.getNextSeq();
    count++;
  }
  return list;
}

ProteoformPtr readFastaToProteoformPtr(FastaIndexReaderPtr reader_ptr,
                                       const std::string &seq_name,
                                       const std::string &seq_desc,
                                       const ModPtrVec &fix_mod_list) {
  FastaSeqPtr seq_ptr = reader_ptr->readFastaSeq(seq_name, seq_desc);
  if (seq_ptr != nullptr) {
    return geneDbProteoformPtr(seq_ptr, fix_mod_list);
  } else {
    return ProteoformPtr(nullptr);
  }
}

}  // namespace toppiceoform_factory

}  // namespace toppic

