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

#include <cmath>
#include <vector>
#include <string>

#include "base/residue_util.hpp"
#include "tdgf/tdgf_util.hpp"

namespace prot {

namespace tdgf_util {

void updateResidueCounts(const ResiduePtrVec &residue_list,
                         std::vector<double> &counts,
                         ProteoformPtr prot_ptr) {
  ResSeqPtr seq_ptr = prot_ptr->getResSeqPtr();
  for (int i = 0; i < seq_ptr->getLen(); i++) {
    ResiduePtr res_ptr = seq_ptr->getResiduePtr(i);
    int pos = residue_util::findResidue(residue_list, res_ptr);
    if (pos >= 0) {
      // found
      counts[pos] = counts[pos]+1;
    }
  }
}

ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                              const std::vector<double> &counts) {
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

void updateNTermResidueCounts(ResiduePtrVec &residue_list, std::vector<double> &counts,
                              const ProteoformPtrVec &mod_proteo_ptrs) {
  for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
    ResSeqPtr seq_ptr = mod_proteo_ptrs[i]->getResSeqPtr();
    if (seq_ptr->getLen() >= 1) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
      int pos = residue_util::findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found
        counts[pos] = counts[pos]+1;
      } else {
        residue_list.push_back(res_ptr);
        counts.push_back(1);
      }
    }
  }
}

int computeAvgLength(const ResFreqPtrVec &residue_ptrs, double convert_ratio) {
  double mass_sum = 0;
  double freq_sum = 0;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    freq_sum = freq_sum + residue_ptrs[i]->getFreq();
    mass_sum = mass_sum + residue_ptrs[i]->getFreq() * residue_ptrs[i]->getMass();
  }
  return static_cast<int>(std::round(mass_sum / freq_sum * convert_ratio));
}

} // namespace tdgf_util

}  // namespace prot
