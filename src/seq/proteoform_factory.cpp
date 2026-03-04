//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#include "common/util/logger.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/residue_util.hpp"
#include "common/base/trunc_util.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/prot_mod_util.hpp"
#include "seq/residue_seq.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace proteoform_factory {

ProteoformPtr geneDbProteoformPtr(FastaSeqPtr fasta_seq_ptr, ModPtrVec fix_mod_list, 
                                  int start_pos) {
  if (fasta_seq_ptr == nullptr) {
    return ProteoformPtr(nullptr);
  }
  ProtModPtr none_prot_mod_ptr = ProtModBase::getProtModPtr_NONE();
  ResiduePtrVec residue_ptrs 
      = residue_util::convertStrToResiduePtrVec(fasta_seq_ptr->getAcidPtmPairVec());
  int end_pos = static_cast<int>(residue_ptrs.size()) - 1;

  MassShiftPtrVec shift_list;
  // add input ptms;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    if (residue_ptrs[i]->getPtmPtr() != PtmBase::getEmptyPtmPtr()) {
      ResiduePtr ori_residue 
          = ResidueBase::getBaseResiduePtr(residue_ptrs[i]->getAminoAcidPtr());
      ModPtr mod_ptr = ModBase::getBaseModPtr(ori_residue, residue_ptrs[i], ModType::SIDE_CHAIN);
      AlterPtr alter_ptr = std::make_shared<Alter>(i, i + 1, AlterType::INPUT, 
                                                   mod_ptr->getShift(), mod_ptr);
      MassShiftPtr shift_ptr = std::make_shared<MassShift>(alter_ptr);
      shift_list.push_back(shift_ptr);
    }
  }

  // add fixed ptms;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    for (size_t j = 0; j < fix_mod_list.size(); j++) {
      if (residue_ptrs[i] == fix_mod_list[j]->getOriResiduePtr()) {
        residue_ptrs[i] = fix_mod_list[j]->getModResiduePtr();
        AlterPtr alter_ptr = std::make_shared<Alter>(i, i + 1, AlterType::FIXED, 
                                                     fix_mod_list[j]->getShift(), 
                                                     fix_mod_list[j]);
        MassShiftPtr shift_ptr = std::make_shared<MassShift>(alter_ptr);
        shift_list.push_back(shift_ptr);
        break;
      }
    }
  }

  ResSeqPtr res_seq_ptr = std::make_shared<ResidueSeq>(residue_ptrs);
  return std::make_shared<Proteoform>(fasta_seq_ptr, none_prot_mod_ptr, start_pos,
                                      end_pos, res_seq_ptr, shift_list);
}

ProteoformPtr geneDbProteoformPtr(FastaSeqPtr fasta_seq_ptr, ModPtrVec fix_mod_list) {
  return geneDbProteoformPtr(fasta_seq_ptr, fix_mod_list, 0);
}

ProteoformPtr geneProtModProteoform(ProteoformPtr db_form_ptr, ProtModPtr prot_mod_ptr) {
  // check if the proteoform can be truncated
  ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
  bool valid_mod = prot_mod_util::allowMod(prot_mod_ptr, db_res_seq_ptr->getResidues());
  if (!valid_mod) {
    //LOG_DEBUG("NO valid mod");
    return ProteoformPtr(nullptr);
  }

  TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
  // first residue index
  int start = trunc_ptr->getTruncLen();
  // last residue index
  int end = db_form_ptr->getLen() - 1;
  // copy input changes
  MassShiftPtrVec ori_shift_ptrs = db_form_ptr->getMassShiftPtrVec();
  MassShiftPtrVec shift_ptrs;
  for (size_t i = 0; i < ori_shift_ptrs.size(); i++) {
    if (ori_shift_ptrs[i]->getLeftBpPos() >= start
        && ori_shift_ptrs[i]->getRightBpPos() <= end + 1) {
      MassShiftPtr shift_ptr = std::make_shared<MassShift>(ori_shift_ptrs[i], start);
      shift_ptrs.push_back(shift_ptr);
    } 
  }

  ResiduePtrVec ori_vec = db_res_seq_ptr->getResidues();
  ResiduePtrVec new_vec(ori_vec.begin()+start, ori_vec.begin()+end+1);
  // apply mod
  ModPtr n_term_mod_ptr = prot_mod_ptr->getModPtr();
  if (!ModBase::isNTermNoneModPtr(n_term_mod_ptr)) {
    int mod_pos = prot_mod_ptr->getModPos();
    int new_pos = mod_pos - start;
    AlterPtr alter_ptr = std::make_shared<Alter>(new_pos, new_pos + 1, 
                                                 AlterType::PROTEIN_VARIABLE,
                                                 n_term_mod_ptr->getShift(), 
                                                 n_term_mod_ptr);
    MassShiftPtr shift_ptr = std::make_shared<MassShift>(alter_ptr); 
    shift_ptrs.push_back(shift_ptr);
  }
  ResSeqPtr seq_ptr = std::make_shared<ResidueSeq>(new_vec, n_term_mod_ptr);

  FastaSeqPtr fasta_seq_ptr = db_form_ptr->getFastaSeqPtr();

  return std::make_shared<Proteoform>(db_form_ptr->getFastaSeqPtr(), prot_mod_ptr, start,
                                      end, seq_ptr, shift_ptrs);
}

ProteoformPtrVec geneProtModProteoform(ProteoformPtr proteo_ptr, 
                                       const ProtModPtrVec &prot_mods) {
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

ProteoformPtr geneSubProteoform(ProteoformPtr proteoform_ptr, 
                                FastaSeqPtr fasta_seq_ptr,
                                int residue_start, int residue_end) {

  ResiduePtrVec ori_residues = proteoform_ptr->getResSeqPtr()->getResidues();
  ResiduePtrVec residues(ori_residues.begin() + residue_start, 
                         ori_residues.begin() + residue_end + 1);
  // if the sub-proteoform starts from the first residue, keep the N-terminal mod  
  // otherwise, set it to none
  ModPtr n_term_mod_ptr = ModBase::getNTermNoneModPtr();
  if (residue_start == 0) {
    n_term_mod_ptr = proteoform_ptr->getResSeqPtr()->getNModPtr();
  }
  ResSeqPtr seq_ptr = std::make_shared<ResidueSeq>(residues, n_term_mod_ptr);

  MassShiftPtrVec ori_shift_list = proteoform_ptr->getMassShiftPtrVec();
  MassShiftPtrVec shift_list;

  for (size_t i = 0; i < ori_shift_list.size(); i++) {
    if (ori_shift_list[i]->getLeftBpPos() >= residue_start
        && ori_shift_list[i]->getRightBpPos() <= residue_end + 1) {
      MassShiftPtr shift_ptr = std::make_shared<MassShift>(ori_shift_list[i], residue_start); 
      shift_list.push_back(shift_ptr);
    } 
  }

  ProtModPtr prot_mod_ptr = proteoform_ptr->getProtModPtr();
  if (residue_start + proteoform_ptr->getStartPos() > 1) {
    prot_mod_ptr = ProtModBase::getProtModPtr_NONE();
  }

  return std::make_shared<Proteoform>(fasta_seq_ptr, prot_mod_ptr,
                                      residue_start + proteoform_ptr->getStartPos(),
                                      residue_end + proteoform_ptr->getStartPos(),
                                      seq_ptr, shift_list);
}

ProteoformPtrVec readFastaToProteoformPtrVec(const std::string &file_name,
                                             const ModPtrVec &fix_mod_list) {
  FastaReader reader(file_name);
  ProteoformPtrVec list;
  FastaSeqPtr seq_ptr = reader.getNextSeq();
  int count = 0;
  while (seq_ptr != nullptr) {
    ProteoformPtr proteo_ptr = geneDbProteoformPtr(seq_ptr, fix_mod_list);
    list.push_back(proteo_ptr);
    seq_ptr = reader.getNextSeq();
    count++;
  }
  LOG_DEBUG("Seq num: " << count);
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

