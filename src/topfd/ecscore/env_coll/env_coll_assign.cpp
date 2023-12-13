//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include <numeric>
#include <algorithm>

#include "common/util/logger.hpp"

#include "ms/spec/peak_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/msmap/ms_map_row_header.hpp"

#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/env_coll/env_coll_assign.hpp"

namespace toppic {

namespace env_coll_assign {

bool checkEnvColl(MsHeaderPtr header_ptr, EnvCollPtrVec &env_coll_list) {
  int ms1_id = header_ptr->getMsOneId();
  double prec_win_bgn = header_ptr->getPrecWinBegin();
  double prec_win_end = header_ptr->getPrecWinEnd();

  SpecFeaturePtrVec new_spec_feats;
  for (size_t coll_id = 0; coll_id < env_coll_list.size(); coll_id++) {
    // check retention time range
    if (ms1_id < env_coll_list[coll_id]->getStartSpecId() 
        || ms1_id > env_coll_list[coll_id]->getEndSpecId()) {
      continue;
    }
    EnvSetPtrVec env_sets = env_coll_list[coll_id]->getEnvSetList();
    double feature_mono_mass = env_coll_list[coll_id]->getMonoNeutralMass();
    double feature_avg_mass = EnvBase::convertMonoMassToAvgMass(feature_mono_mass); 
    for (size_t env_set_id = 0; env_set_id < env_sets.size(); env_set_id++) {
      EnvSetPtr env_set_ptr = env_sets[env_set_id];
      double mz = peak_util::compMz(feature_avg_mass, env_set_ptr->getCharge());
      // check precsor window
      // this part needs to be improved to consider only peaks in the window
      if (mz < prec_win_bgn || mz > prec_win_end) {
        continue;
      }
      // check retentime time range
      if (ms1_id < env_set_ptr->getStartSpecId() 
          || ms1_id > env_set_ptr->getEndSpecId()) {
        continue;
      }
      // get intensity information
      std::vector<double> env_intes = env_set_ptr->getXicPtr()->getAllPeakInteList();
      if (env_intes.size() == 0) {
        LOG_WARN("Empty envelope intensity list!");
        continue; 
      }
      return true;
    }
  }
  return false;
}

bool getHighestInteEnvColl(FracFeaturePtrVec &frac_features, EnvCollPtrVec &env_coll_list,
                           MsHeaderPtr header_ptr, double score_thresh, 
                           SpecFeaturePtrVec &ms2_features) {
  int ms1_id = header_ptr->getMsOneId();
  double prec_win_bgn = header_ptr->getPrecWinBegin();
  double prec_win_end = header_ptr->getPrecWinEnd();

  SpecFeaturePtrVec new_spec_feats;
  for (size_t coll_id = 0; coll_id < env_coll_list.size(); coll_id++) {
    FracFeaturePtr frac_feature_ptr = frac_features[coll_id];
    // check threshold
    if (frac_feature_ptr->getEcScore() < score_thresh) {
      continue;
    }
    // check retention time range
    if (ms1_id < env_coll_list[coll_id]->getStartSpecId() 
        || ms1_id > env_coll_list[coll_id]->getEndSpecId()) {
      continue;
    }
    EnvSetPtrVec env_sets = env_coll_list[coll_id]->getEnvSetList();
    double feature_mono_mass = env_coll_list[coll_id]->getMonoNeutralMass();
    double feature_avg_mass = EnvBase::convertMonoMassToAvgMass(feature_mono_mass); 
    for (size_t env_set_id = 0; env_set_id < env_sets.size(); env_set_id++) {
      EnvSetPtr env_set_ptr = env_sets[env_set_id];
      double mz = peak_util::compMz(feature_avg_mass, env_set_ptr->getCharge());
      // check precsor window
      // this part needs to be improved to consider only peaks in the window
      if (mz < prec_win_bgn || mz > prec_win_end) {
        continue;
      }
      // check retentime time range
      if (ms1_id < env_set_ptr->getStartSpecId() 
          || ms1_id > env_set_ptr->getEndSpecId()) {
        continue;
      }
      // get intensity information
      int inte_idx = ms1_id - env_set_ptr->getStartSpecId();
      std::vector<double> env_intes = env_set_ptr->getXicPtr()->getAllPeakInteList();
      if (env_intes.size() == 0) {
        LOG_WARN("Empty envelope intensity list!");
        continue; 
      }

      int prec_charge = env_set_ptr->getCharge();
      double prec_mono_mz = peak_util::compMz(feature_mono_mass, prec_charge);
      double prec_avg_mz = peak_util::compMz(feature_avg_mass, prec_charge); 
      double prec_inte = env_intes[inte_idx];
      if (prec_inte < 0) {
        prec_inte = 0;
      }
      SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(header_ptr,
                                                                 frac_feature_ptr, 
                                                                 prec_mono_mz, 
                                                                 prec_avg_mz, 
                                                                 prec_charge,
                                                                 prec_inte);
      new_spec_feats.push_back(ms2_feature);
      frac_feature_ptr->setHasMs2Spec(true);
      break;
    }
  }
  if (new_spec_feats.size() > 0) {
    // sort by intensity
    std::sort(new_spec_feats.begin(), new_spec_feats.end(),
              SpecFeature::cmpPrecInteDec);
    ms2_features.insert(std::end(ms2_features), std::begin(new_spec_feats), 
                        std::end(new_spec_feats)); 
    return true;
  }
  else {
    return false;
  }
}

void assignEnvColls(FracFeaturePtrVec &frac_feature_list,
                    EnvCollPtrVec &env_coll_list,
                    MsHeaderPtr2D &ms2_header_ptr_2d,
                    SpecFeaturePtrVec &ms2_feature_list, 
                    double score_cutoff) {
  for (size_t ms1_idx = 0; ms1_idx < ms2_header_ptr_2d.size(); ms1_idx++) {
    for (size_t i = 0; i < ms2_header_ptr_2d[ms1_idx].size(); i++) {
      MsHeaderPtr header_ptr = ms2_header_ptr_2d[ms1_idx][i];
      bool assigned = getHighestInteEnvColl(frac_feature_list, env_coll_list,  
                                            header_ptr, score_cutoff, ms2_feature_list);
      if (!assigned) {
        // lower the cutoff to 0. This step is necessary when low intensity
        // features are added.
        score_cutoff = 0;
        assigned = getHighestInteEnvColl(frac_feature_list, env_coll_list,  
                                         header_ptr, score_cutoff, ms2_feature_list);
      }
      if (!assigned) {
        LOG_INFO("Scan " << header_ptr->getFirstScanNum() << " does not have MS1 feature!");
      }
    }
  }
  std::sort(ms2_feature_list.begin(), ms2_feature_list.end(), SpecFeature::cmpSpecIdInc);
}

}

}
