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

#include <vector>
#include <algorithm>
#include <cmath>

#include "common/base/residue_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/prot_mod_util.hpp"
#include "seq/proteoform_util.hpp"

namespace toppic {
namespace proteoform_util {

ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts;
  ResiduePtrVec residue_list;
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
    if (seq_ptr->getLen() >= 1) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
      int pos = residue_util::findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found
        counts[pos] = counts[pos] + 1;
      } else {
        residue_list.push_back(res_ptr);
        counts.push_back(1);
      }
    }
  }

  double sum = 0;
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr = std::make_shared<ResidueFreq>(residue_list[i]->getAminoAcidPtr(),
                                                            residue_list[i]->getPtmPtr(),
                                                            counts[i] / sum);
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                              const ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts(residue_list.size(), 0.0);
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
    for (int j = 0; j < seq_ptr->getLen(); j++) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(j);
      int pos = residue_util::findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found
        counts[pos] = counts[pos]+1;
      }
    }
  }

  double sum = 0;
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr = std::make_shared<ResidueFreq>(residue_list[i]->getAminoAcidPtr(),
                                                            residue_list[i]->getPtmPtr(),
                                                            counts[i] / sum);
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

bool isSameSeqAndMass(ProteoformPtr a, ProteoformPtr b, double ppo) {
  if (a->getSeqName() != b->getSeqName()) {
    return false;
  }

  if (a->getStartPos() != b->getStartPos()) {
    return false;
  }

  if (a->getEndPos() != b->getEndPos()) {
    return false;
  }

  double thresh = a->getMass() * ppo;

  if (std::abs(a->getMass() -b->getMass())> thresh) {
    return false;
  }

  return true;
}

bool isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b, double ppo) {
  if (!isSameSeqAndMass(a, b, ppo)) {
    return false;
  }

  if (a->getMassShiftNum() != b->getMassShiftNum()) {
    return false;
  }

  double shift_tolerance = a->getResSeqPtr()->getSeqMass() * ppo;
  // sort mass shifts
  MassShiftPtrVec a_shift_vec = a->getMassShiftPtrVec();
  MassShiftPtrVec b_shift_vec = b->getMassShiftPtrVec();
  std::sort(a_shift_vec.begin(), a_shift_vec.end(), MassShift::cmpPosInc);
  std::sort(b_shift_vec.begin(), b_shift_vec.end(), MassShift::cmpPosInc);
  for (int i = 0; i < a->getMassShiftNum(); i++) {
    MassShiftPtr ac = a_shift_vec[i];
    MassShiftPtr bc = b_shift_vec[i];
    if (ac->getRightBpPos() <= bc->getLeftBpPos() || bc->getRightBpPos() <= ac->getLeftBpPos()) {
      return false;
    }
    if (std::abs(ac->getMassShift()-bc->getMassShift()) > shift_tolerance) {
      return false;
    }
  }
  return true;
}

ProteoformPtrVec2D divideProteoIntoBlocks(const ProteoformPtrVec &proteo_ptrs, int db_block_size) {
  size_t start_idx = 0;
  size_t proteo_idx = 0;
  int block_len = 0;
  ProteoformPtrVec2D proteo_blocks;
  while (proteo_idx < proteo_ptrs.size()) {
    int proteo_len = proteo_ptrs[proteo_idx]->getResSeqPtr()->getLen();
    if (block_len + proteo_len < db_block_size) {
      block_len = block_len + proteo_len;
      proteo_idx++;
    } else {
      size_t end_idx = proteo_idx;
      ProteoformPtrVec proteo_in_block;
      for (size_t i = start_idx; i <= end_idx; i++) {
        proteo_in_block.push_back(proteo_ptrs[i]);
      }
      proteo_blocks.push_back(proteo_in_block);
      start_idx = end_idx +1;
      proteo_idx = end_idx +1;
      block_len = 0;
    }
  }
  /* last block */
  if (start_idx < proteo_ptrs.size()) {
    ProteoformPtrVec proteo_in_block;
    for (size_t i = start_idx; i < proteo_ptrs.size(); i++) {
      proteo_in_block.push_back(proteo_ptrs[i]);
    }
    proteo_blocks.push_back(proteo_in_block);
  }
  return proteo_blocks;
}

std::vector<double> getNTermShift(ProteoformPtr db_form_ptr,
                                  const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<double> shifts;
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
    bool valid = prot_mod_util::allowMod(prot_mod_ptrs[i], db_res_seq_ptr->getResidues());
    // LOG_DEBUG("valid " << valid << " shift " << prot_mod_ptrs[i]->getProtShift());
    if (valid) {
      shifts.push_back(prot_mod_ptrs[i]->getProtShift());
    }
  }
  return shifts;
}

std::vector<double> getNTermAcets(ProteoformPtr db_form_ptr,
                                  const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<double> shifts;
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    // check if it is acetylation
    if (prot_mod_ptrs[i]->getModPtr()->getModResiduePtr()->getPtmPtr()
        == PtmBase::getPtmPtr_Acetylation()) {
      ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
      bool valid = prot_mod_util::allowMod(prot_mod_ptrs[i], db_res_seq_ptr->getResidues());
      if (valid) {
        shifts.push_back(prot_mod_ptrs[i]->getProtShift());
      }
    }
  }
  return shifts;
}

std::vector<std::vector<double> > getNTermShift2D(const ProteoformPtrVec & db_form_ptr_vec,
                                                  const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<std::vector<double> > shifts_2d;
  for (size_t i = 0; i < db_form_ptr_vec.size(); i++) {
    std::vector<double> shifts = getNTermShift(db_form_ptr_vec[i], prot_mod_ptrs);
    shifts_2d.push_back(shifts);
  }
  return shifts_2d;
}

std::vector<std::vector<double> > getNTermAcet2D(const ProteoformPtrVec & db_form_ptr_vec,
                                                 const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<std::vector<double> > shifts_2d;
  for (size_t i = 0; i < db_form_ptr_vec.size(); i++) {
    std::vector<double> shifts = getNTermAcets(db_form_ptr_vec[i], prot_mod_ptrs);
    shifts_2d.push_back(shifts);
  }
  return shifts_2d;
}

ProteoformPtr geneDbProteoformPtr(FastaSubSeqPtr fasta_seq_ptr, ModPtrVec fix_mod_list, int start_pos) {
  if (fasta_seq_ptr == nullptr) {
    return nullptr;
  }
  ProtModPtr none_prot_mod_ptr = ProtModBase::getProtModPtr_NONE();
  ResiduePtrVec residue_ptrs = residue_util::convertStrToResiduePtrVec(fasta_seq_ptr->getAcidPtmPairVec());
  int end_pos = start_pos + static_cast<int>(residue_ptrs.size()) - 1;

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

}  // namespace toppiceoform_util
}  // namespace toppic

