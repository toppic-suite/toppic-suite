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

  PrsmPtr processOnePtm(PrsmPtr prsm);

  PrsmPtr processTwoPtm(PrsmPtr prsm);

  ProteoformPtr processOneKnownPtm(PrsmPtr prsm);

  ProteoformPtr processTwoKnownPtm(PrsmPtr prsm);

  bool modifiable(ProteoformPtr proteoform_ptr, int i, PtmPtr ptm_ptr);

  void compOnePtmScr(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                     std::vector<double> &scr_vec, double & raw_scr, PtmPtrVec & ptm_vec);

  void compTwoPtmScr(ProteoformPtr proteoform, int num_match,
                     const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                     double & raw_scr, PtmPairVec & ptm_pair_vec);

  double dpTwoPtmScr(ProteoformPtr proteoform, int h, const ExtendMsPtrVec & extend_ms_ptr_vec,
                     double prec_mass, double mass1, double mass2, PtmPtr ptm1, PtmPtr ptm2);

  ProteoformPtr onePtmTermAdjust(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                                 double & shift_mass, double err);

  ProteoformPtr twoPtmTermAdjust(ProteoformPtr proteoform, int num_match,
                                 const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                                 double & mass1, double & mass2);

  ProteoformPtr compSplitPoint(ProteoformPtr proteoform, int h, const ExtendMsPtrVec & extend_ms_ptr_vec,
                               double prec_mass);

  void getNtermTruncRange(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                          int & min, int & max);

  void getCtermTruncRange(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                          int & min, int & max);

  LocalMngPtr mng_ptr_;

  // Single PTMs
  PtmPtrVec ptm_vec_;
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
