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


#include <algorithm>
#include <cmath>

#include "base/residue_util.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"
#include "base/proteoform_util.hpp"

namespace prot {

ResFreqPtrVec ProteoformUtil::compNTermResidueFreq(
    const ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts;
  ResiduePtrVec residue_list;
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
    if (seq_ptr->getLen() >= 1) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
      int pos = ResidueUtil::findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found
        counts[pos] = counts[pos]+1;
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
    ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(),
                                            residue_list[i]->getPtmPtr(),
                                            counts[i]/sum));
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}


ResFreqPtrVec ProteoformUtil::compResidueFreq(const ResiduePtrVec &residue_list,
                                              const ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts(residue_list.size(), 0.0);
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
    for (int j = 0; j < seq_ptr->getLen(); j++) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(j);
      int pos = ResidueUtil::findResidue(residue_list, res_ptr);
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
    ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(),
                                            residue_list[i]->getPtmPtr(),
                                            counts[i]/sum));
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

bool ProteoformUtil::isSameSeqAndMass(ProteoformPtr a, ProteoformPtr b, 
                                      double ppo) {
  if(a->getSeqName() != b->getSeqName()) {
    return false;
  }
  if(a->getStartPos() != b->getStartPos()) {
    return false;
  }
  if(a->getEndPos() != b->getEndPos()) {
    return false;
  }
  double thresh = a->getMass() * ppo;
  if(std::abs(a->getMass() -b->getMass())> thresh) {
    return false;
  }
  return true;
}

bool ProteoformUtil::isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b,
                                                   double ppo) {
  if(!isSameSeqAndMass(a,b,ppo)) {
    return false;
  }
  if(a->getChangePtrVec().size() != b->getChangePtrVec().size()) {
    return false;
  }
  double shift_tolerance = a->getResSeqPtr()->getSeqMass()*ppo;
  // sort changes
  ChangePtrVec a_change_vec = a->getChangePtrVec();
  ChangePtrVec b_change_vec = b->getChangePtrVec();
  std::sort(a_change_vec.begin(),a_change_vec.end(),Change::cmpPosInc);
  std::sort(b_change_vec.begin(),b_change_vec.end(),Change::cmpPosInc);
  for(size_t i=0; i< a->getChangePtrVec().size(); i++) {
    ChangePtr ac = a_change_vec[i];
    ChangePtr bc = b_change_vec[i];
    if(ac->getRightBpPos() <= bc->getLeftBpPos() || bc->getRightBpPos() <= ac->getLeftBpPos()) {
      return false;
    }
    if(std::abs(ac->getMassShift()-bc->getMassShift()) > shift_tolerance) {
      return false;
    }
  }
  return true;
}

ProteoformPtrVec2D ProteoformUtil::divideProteoIntoBlocks(
    const ProteoformPtrVec &proteo_ptrs, int db_block_size) {
  size_t start_idx = 0;
  size_t proteo_idx =0;
  int block_len=0;
  ProteoformPtrVec2D proteo_blocks;
  while(proteo_idx < proteo_ptrs.size()) {
    int proteo_len = proteo_ptrs[proteo_idx]->getResSeqPtr()->getLen();
    if(block_len + proteo_len < db_block_size) {
      block_len = block_len + proteo_len;
      proteo_idx++;
    } else {
      size_t end_idx = proteo_idx;
      ProteoformPtrVec proteo_in_block;
      for(size_t i = start_idx; i<=end_idx; i++) {
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
    for(size_t i = start_idx; i < proteo_ptrs.size(); i++) {
      proteo_in_block.push_back(proteo_ptrs[i]);
    }
    proteo_blocks.push_back(proteo_in_block);
  }
  return proteo_blocks;
}

std::vector<double> ProteoformUtil::getNTermShift(ProteoformPtr db_form_ptr,
                                               const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<double> shifts;
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
    bool valid = ProtModUtil::allowMod(prot_mod_ptrs[i], 
                                       db_res_seq_ptr->getResidues());
    //LOG_DEBUG("valid " << valid << " shift " << prot_mod_ptrs[i]->getProtShift());
    if (valid) {
      shifts.push_back(prot_mod_ptrs[i]->getProtShift());
    }
  }
  return shifts;
}

std::vector<double> ProteoformUtil::getNTermAcets(ProteoformPtr db_form_ptr,
                                               const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<double> shifts;
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    // check if it is acetylation
    if (prot_mod_ptrs[i]->getModPtr()->getModResiduePtr()->getPtmPtr() 
        == PtmBase::getPtmPtr_Acetylation()) {
      ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
      bool valid = ProtModUtil::allowMod(prot_mod_ptrs[i], 
                                         db_res_seq_ptr->getResidues());
      if (valid) {
        shifts.push_back(prot_mod_ptrs[i]->getProtShift());
      }
    }
  }
  return shifts;
}

std::vector<std::vector<double>> ProteoformUtil::getNTermShift2D(
    ProteoformPtrVec db_form_ptr_vec, const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<std::vector<double>> shifts_2d;
  for (size_t i = 0; i < db_form_ptr_vec.size(); i++) {
    std::vector<double> shifts = getNTermShift(db_form_ptr_vec[i], prot_mod_ptrs);
    shifts_2d.push_back(shifts);
  }
  return shifts_2d;
}

std::vector<std::vector<double>> ProteoformUtil::getNTermAcet2D(
    ProteoformPtrVec db_form_ptr_vec, const ProtModPtrVec &prot_mod_ptrs) {
  std::vector<std::vector<double>> shifts_2d;
  for (size_t i = 0; i < db_form_ptr_vec.size(); i++) {
    std::vector<double> shifts = getNTermAcets(db_form_ptr_vec[i], prot_mod_ptrs);
    shifts_2d.push_back(shifts);
  }
  return shifts_2d;
}

} /* namespace prot */

