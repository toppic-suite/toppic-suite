//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "merge/feature_prsm.hpp"

namespace toppic {

FeaturePrsm::FeaturePrsm(PrsmStrPtr prsm) {
  prot_name_ = prsm->getSeqName();
  prot_desc_ = prsm->getSeqDesc();
  first_residue_ = prsm->getProteoformStartPos();
  last_residue_ = prsm->getProteoformEndPos();
  proteoform_ = prsm->getProteoformMatchSeq();
  ms2_id_ = prsm->getSpectrumId();
  prec_mass_ = prsm->getOriPrecMass();
  proteo_id_ = prsm->getProteoClusterId();
  proteo_inte_ = prsm->getProteoInte();
  min_time_ = prsm->getFracFeatureMinTime();
  max_time_ = prsm->getFracFeatureMaxTime();
  apex_time_ = prsm->getFracFeatureApexTime();
}

}

