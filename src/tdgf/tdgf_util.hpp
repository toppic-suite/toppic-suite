//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_TDGF_UTIL_HPP_
#define PROT_TDGF_UTIL_HPP_

#include "base/proteoform.hpp"
#include "base/residue_freq.hpp"
#include "tdgf/tdgf_mng.hpp"

namespace prot {

class TdgfUtil {
 public:
  static void updateNTermResidueCounts(ResiduePtrVec &residue_list, 
                                       std::vector<double> &counts,
                                       const ProteoformPtrVec &mod_proteo_ptrs);


  static void updateResidueCounts(const ResiduePtrVec &residue_list, 
                                  std::vector<double> &counts,
                                  ProteoformPtr prot_ptr);

  static ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list, 
                                       const std::vector<double> &counts);

  static int computeAvgLength(const ResFreqPtrVec &residue_ptrs, double convert_ratio);
};

}

#endif
