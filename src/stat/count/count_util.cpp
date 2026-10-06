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

#include "stat/count/count_util.hpp"

#include <cmath>

#include "common/base/residue_util.hpp"

namespace toppic {

namespace count_util {

// A residue not in the list, one with a PTM written in the FASTA sequence
// (M[Oxidation]), is added to it, as the residue of a fixed modification is
// in it; it used to be left out of the frequencies.
void updateResidueCounts(ResiduePtrVec& residue_list,
                         std::vector<double>& counts,
                         const ProteoformPtr& prot_ptr) {
  ResSeqPtr seq_ptr = prot_ptr->getResSeqPtr();
  for (int i = 0; i < seq_ptr->getLen(); i++) {
    ResiduePtr res_ptr = seq_ptr->getResiduePtr(i);
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

ResFreqPtrVec compResidueFreq(const ResiduePtrVec& residue_list,
                              const std::vector<double>& counts) {
  double sum = 0;
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr = std::make_shared<ResidueFreq>(
        residue_list[i]->getAminoAcidPtr(), residue_list[i]->getPtmPtr(),
        counts[i] / sum);
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

void updateNTermResidueCounts(ResiduePtrVec& residue_list,
                              std::vector<double>& counts,
                              const ProteoformPtrVec& mod_proteo_ptrs) {
  for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
    ResSeqPtr seq_ptr = mod_proteo_ptrs[i]->getResSeqPtr();
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
}

int computeAvgLength(const ResFreqPtrVec& residue_ptrs, double convert_ratio) {
  double mass_sum = 0;
  double freq_sum = 0;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    freq_sum = freq_sum + residue_ptrs[i]->getFreq();
    mass_sum =
        mass_sum + residue_ptrs[i]->getFreq() * residue_ptrs[i]->getMass();
  }
  return static_cast<int>(std::round(mass_sum / freq_sum * convert_ratio));
}

ResFreqPtrVec geneModResFreqs(const ResFreqPtrVec& residue_ptrs,
                              const ModPtrVec& mod_ptr_vec) {
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    // the modified residues of this residue, each once: a residue is
    // modified as the variable PTM search modifies it, when it is the
    // modification's original residue (the same amino acid without a PTM),
    // so the residue of a fixed modification (C57) is not
    ResiduePtrVec mod_res_ptrs;
    for (size_t j = 0; j < mod_ptr_vec.size(); j++) {
      const ResiduePtr& mod_res_ptr = mod_ptr_vec[j]->getModResiduePtr();
      if (mod_ptr_vec[j]->getOriResiduePtr()->isSame(residue_ptrs[i]) &&
          residue_util::findResidue(mod_res_ptrs, mod_res_ptr) < 0) {
        mod_res_ptrs.push_back(mod_res_ptr);
      }
    }
    double freq = residue_ptrs[i]->getFreq() / (mod_res_ptrs.size() + 1);
    res_freq_list.push_back(std::make_shared<ResidueFreq>(
        residue_ptrs[i]->getAminoAcidPtr(), residue_ptrs[i]->getPtmPtr(),
        freq));
    for (size_t j = 0; j < mod_res_ptrs.size(); j++) {
      res_freq_list.push_back(std::make_shared<ResidueFreq>(
          mod_res_ptrs[j]->getAminoAcidPtr(), mod_res_ptrs[j]->getPtmPtr(),
          freq));
    }
  }
  return res_freq_list;
}

}  // namespace count_util

}  // namespace toppic
