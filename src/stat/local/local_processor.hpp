//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_STAT_LOCAL_PROCESSOR_HPP_
#define TOPPIC_STAT_LOCAL_PROCESSOR_HPP_

#include "common/base/ptm.hpp"
#include "prsm/prsm.hpp"
#include "stat/local/local_mng.hpp"

namespace toppic {

class LocalProcessor {
 public:
  explicit LocalProcessor(LocalMngPtr mng_ptr);

  void process();

 private:
  void init();

  PrsmPtr processOnePrsm(PrsmPtr prsm);

  PrsmPtr processOneMassShift(PrsmPtr prsm);

  PrsmPtr processTwoMassShifts(PrsmPtr prsm);

  ProteoformPtr processOneKnownPtm(PrsmPtr prsm);

  ProteoformPtr processTwoKnownPtms(PrsmPtr prsm);

  int compOnePtmScr(ProteoformPtr base_form_ptr, 
                    const ExtendMsPtrVec & extend_ms_ptr_vec,
                    PtmPtr ptm_ptr); 

  ProteoformPtr onePtmLocalize(ProteoformPtr base_form_ptr, 
                               const ExtendMsPtrVec & extend_ms_ptr_vec,
                               double shift_mass, int match_score,
                               PtmPtr ptm_ptr);

  int compTwoPtmScr(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec, 
                    PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2); 

  ProteoformPtr twoPtmLocalize(ProteoformPtr form_ptr, const ExtendMsPtrVec &extend_ms_ptr_vec,
                               int match_score, PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2); 

  bool modifiable(ProteoformPtr proteoform_ptr, int i, PtmPtr ptm_ptr);

  LocalMngPtr mng_ptr_;

  // Single PTMs
  PtmPtrVec ptm_ptr_vec_;
  // Ptm Pairs
  PtmPairVec ptm_pair_vec_;
  // N-terminal modification list
  ModPtrVec mod_list_N_;
  // C-terminal modification list
  ModPtrVec mod_list_C_;
  // Modification list for any positions
  ModPtrVec mod_list_any_;
};

typedef std::shared_ptr<LocalProcessor> LocalProcessorPtr;

}  // namespace toppic

#endif
