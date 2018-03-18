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


#ifndef PROT_LOCAL_PROCESSOR_HPP_
#define PROT_LOCAL_PROCESSOR_HPP_

#include <vector>

#include "base/ptm.hpp"
#include "base/local_anno.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "spec/theo_peak.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "local_mng.hpp"

namespace prot {

class LocalProcessor {
 public:
  explicit LocalProcessor(LocalMngPtr mng_ptr):
      mng_ptr_(mng_ptr),
      ppo_(mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo()),
      theta_(mng_ptr->theta_),
      threshold_(mng_ptr->threshold_),
      beta_(mng_ptr->beta_),
      min_mass_(mng_ptr->min_mass_),
      p1_(mng_ptr->p1_),
      p2_(mng_ptr->p2_) {
        init();
      }

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

  double ppo_;

  double theta_;  // the weight for known/unknown ptm

  double threshold_;  // threshold for MIScore

  double beta_;  // the weight for one/two ptm

  double min_mass_;

  double p1_, p2_;

  PtmPtrVec ptm_vec_;

  PtmPairVec ptm_pair_vec_;

  ModPtrVec mod_list_N_;

  ModPtrVec mod_list_C_;

  ModPtrVec mod_list_any_;
};

typedef std::shared_ptr<LocalProcessor> LocalProcessorPtr;

}  // namespace prot

#endif
