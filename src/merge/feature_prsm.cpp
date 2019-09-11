//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/str_util.hpp"
#include "merge/feature_prsm.hpp"

namespace toppic {

FeaturePrsm::FeaturePrsm(std::string line):
    SampleFeature(line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  if (strs.size() > 10) {
    prot_name_ = strs[10];
    prot_desc_ = strs[11];
    first_residue_ = std::stoi(strs[12]) - 1;
    last_residue_ = std::stoi(strs[13]) - 1;
    proteoform_ = strs[14];
    ms2_id_ = std::stoi(strs[15]);
    prec_mass_ = std::stod(strs[16]);
  }
  else {
    ms2_id_ = -1;
  }
}

void FeaturePrsm::addPrsmInfo(PrsmStrPtr prsm) {
  prot_name_ = prsm->getSeqName();
  prot_desc_ = prsm->getSeqDesc();
  first_residue_ = prsm->getProteoformStartPos();
  last_residue_ = prsm->getProteoformEndPos();
  proteoform_ = prsm->getProteinMatchSeq();
  ms2_id_ = prsm->getSpectrumId();
  prec_mass_ = prsm->getOriPrecMass(); 
}

}

