
#include <algorithm>
#include <cmath>

#include "base/proteoform_util.hpp"
#include "base/residue_util.hpp"

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
  double thresh = a->getResSeqPtr()->getSeqMass() * ppo;
  if(std::abs(a->getResSeqPtr()->getSeqMass()
              -b->getResSeqPtr()->getSeqMass())> thresh) {
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
  std::sort(a_change_vec.begin(),a_change_vec.end(),Change::cmpPosIncrease);
  std::sort(b_change_vec.begin(),b_change_vec.end(),Change::cmpPosIncrease);
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


} /* namespace prot */

