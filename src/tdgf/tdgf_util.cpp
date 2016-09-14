// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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



#include "base/residue_util.hpp"
#include "tdgf/tdgf_util.hpp"

namespace prot {

void TdgfUtil::updateResidueCounts(const ResiduePtrVec &residue_list, 
                                   std::vector<double> &counts,
                                   ProteoformPtr prot_ptr) {
  ResSeqPtr seq_ptr = prot_ptr->getResSeqPtr();    
  for (int i = 0; i < seq_ptr->getLen(); i++) {
    ResiduePtr res_ptr = seq_ptr->getResiduePtr(i);
    int pos = ResidueUtil::findResidue(residue_list, res_ptr);
    if (pos >= 0) {
      // found 
      counts[pos] = counts[pos]+1;
    }
  }
}

ResFreqPtrVec TdgfUtil::compResidueFreq(const ResiduePtrVec &residue_list, 
                                        const std::vector<double> &counts) {
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

void TdgfUtil::updateNTermResidueCounts(ResiduePtrVec &residue_list, std::vector<double> &counts,
                                        const ProteoformPtrVec &mod_proteo_ptrs) {
  for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
    ResSeqPtr seq_ptr = mod_proteo_ptrs[i]->getResSeqPtr();    
    if (seq_ptr->getLen() >= 1) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
      int pos = ResidueUtil::findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found 
        counts[pos] = counts[pos]+1;
      }
      else {
        residue_list.push_back(res_ptr);
        counts.push_back(1);
      }
    }
  }
}

int TdgfUtil::computeAvgLength(const ResFreqPtrVec &residue_ptrs, double convert_ratio) {
  double mass_sum = 0;
  double freq_sum = 0;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    freq_sum = freq_sum + residue_ptrs[i]->getFreq();
    mass_sum = mass_sum + residue_ptrs[i]->getFreq() * residue_ptrs[i]->getMass();
  }
  return (int)std::round(mass_sum/freq_sum * convert_ratio);
}

}
