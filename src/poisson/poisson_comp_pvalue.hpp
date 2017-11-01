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


#ifndef PROT_POSSION_COMP_PVALUE_HPP_
#define PROT_POSSION_COMP_PVALUE_HPP_

#include "tdgf/count_test_num.hpp"
#include "poisson/poisson_mng.hpp"

namespace prot {

class PoissonCompPValue {
 public:
  PoissonCompPValue(ProteoformPtrVec &raw_forms, 
                    ProteoformPtrVec &prot_mod_forms,
                    PoissonMngPtr mng_ptr);

  ExtremeValuePtrVec compExtremeValues(ExtendMsPtr extend_ms, 
                                       PrsmPtrVec &prsms, bool strict);
  void setPValueArray(ExtendMsPtr extend_ms_ptr, PrsmPtrVec &prsms);

 private:
  PoissonMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;
  double residue_avg_len_;

  double compRandMatchProb(double prec_mass, bool is_strict);
  double compConditionProb(double rand_match_prob, ExtendMsPtr extend_ms, PrsmPtr prsm);
};

typedef std::shared_ptr<PoissonCompPValue> PoissonCompPValuePtr;

}

#endif
