//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_STAT_LOCAL_LOCAL_PROCESSOR_HPP_
#define TOPPIC_STAT_LOCAL_LOCAL_PROCESSOR_HPP_

#include "common/base/ptm.hpp"
#include "prsm/prsm.hpp"
#include "stat/local/local_mng.hpp"

namespace toppic {

class LocalProcessor {
 public:
  explicit LocalProcessor(const LocalMngPtr &mng_ptr);

  void process();

 private:
  void init();

  PrsmPtr processOnePrsm(const PrsmPtr &prsm);

  PrsmPtr processOneMassShift(const PrsmPtr &prsm);

  PrsmPtr processTwoMassShifts(const PrsmPtr &prsm);

  ProteoformPtr processOneKnownPtm(const PrsmPtr &prsm);

  ProteoformPtr processTwoKnownPtms(const PrsmPtr &prsm);

  int compOnePtmScr(const ProteoformPtr &base_form_ptr, 
                    const ExtendMsPtrVec & extend_ms_ptr_vec,
                    const PtmPtr &ptm_ptr); 

  ProteoformPtr onePtmLocalize(const ProteoformPtr &base_form_ptr, 
                               const ExtendMsPtrVec & extend_ms_ptr_vec,
                               double shift_mass, int match_score,
                               const PtmPtr &ptm_ptr);

  int compTwoPtmScr(const ProteoformPtr &proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec, 
                    const PtmPtr &ptm_ptr_1, const PtmPtr &ptm_ptr_2); 

  ProteoformPtr twoPtmLocalize(const ProteoformPtr &form_ptr, const ExtendMsPtrVec &extend_ms_ptr_vec,
                               int match_score, const PtmPtr &ptm_ptr_1, const PtmPtr &ptm_ptr_2); 

  bool modifiable(const ProteoformPtr &proteoform_ptr, int i, const PtmPtr &ptm_ptr);

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

using LocalProcessorPtr = std::shared_ptr<LocalProcessor>;

}  // namespace toppic

#endif
