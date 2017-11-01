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


#ifndef PROT_POSSION_PROCESSOR_HPP_
#define PROT_POSSION_PROCESSOR_HPP_

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "poisson/poisson_mng.hpp"
#include "poisson/poisson_comp_pvalue.hpp"

namespace prot {

class PoissonProcessor {
 public:
  PoissonProcessor(PoissonMngPtr mng_ptr);
  void init();
  void process();
  void processOneSpectrum(DeconvMsPtr ms_ptr, PrsmWriter &writer);
 private:
    PoissonMngPtr mng_ptr_;
    ProteoformPtrVec proteoforms_;
    PrsmPtrVec prsms_;
    PoissonCompPValuePtr comp_ptr_;
};

typedef std::shared_ptr<PoissonProcessor> PoissonProcessorPtr;

}

#endif
