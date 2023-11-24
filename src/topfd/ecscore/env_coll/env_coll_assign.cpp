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

bool getNewEnvColl(MsHeaderPtr header_ptr, int ms1_idx, MsMapPtr matrix_ptr,
                   EcscoreParaPtr score_para_ptr, ECScorePtrVec &ecscore_list,
                   EnvCollPtrVec &env_coll_list, DeconvMsPtrVec &ms1_ptr_vec,
                   SeedEnvPtr2D &seed_ptr_2d, FracFeaturePtrVec &frac_feature_list,
                   SpecFeaturePtrVec &ms2_feature_list) {
  // set min match envelope to 1 to accept single scan features
  score_para_ptr->min_match_env_ = 1;
  score_para_ptr->min_match_peak_ = 1;
  score_para_ptr->min_seed_match_peak_ = 0;
  double sn_ratio = 0;
  
  double prec_win_begin = header_ptr->getPrecWinBegin();
  double prec_win_end = header_ptr->getPrecWinEnd();
  SeedEnvPtrVec seed_ptr_list = seed_ptr_2d[ms1_idx];
  SeedEnvPtrVec selected_seed_list;
  LOG_DEBUG("ms1 id  " << ms1_idx << " seed number " << seed_ptr_list.size());
  for (size_t i = 0; i < seed_ptr_list.size(); i++) {
    double ref_mz = seed_ptr_list[i]->getReferMz();
    if (ref_mz > prec_win_begin && ref_mz < prec_win_end) {
      selected_seed_list.push_back(seed_ptr_list[i]);
    }
  }
  SeedEnvPtr seed_ptr;
  if (selected_seed_list.size() > 0) {
    // choose the highest intensity one
    std::sort(selected_seed_list.begin(), selected_seed_list.end(),
              SeedEnv::cmpSeedInteDec);
    seed_ptr = selected_seed_list[0];  
    seed_ptr = seed_env_util::relaxProcessSeedEnvPtr(seed_ptr, matrix_ptr,  
                                                     score_para_ptr, sn_ratio); 
  }
  /* Let us rethink how to handle cases in which precursor peaks are missing
  else {
    // Sometimes the seed envelope reference mz is slightly out of the precursor
    // window. In this case, we use default precursor target mz to generate a
    // seed envelope
    int peak_id = -1;
    int charge = 1;
    double inte = 1;
    double mono_mz = header_ptr->getPrecTargetMz();
    double mono_mass = peak_util::compPeakNeutralMass(mono_mz, charge); 
    DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(ms1_id, peak_id,
                                                          mono_mass, inte, charge);
    seed_ptr = std::make_shared<SeedEnv>(peak_ptr);
  }
  bool valid = seed_env_util::simplePreprocessEnv(matrix_ptr, seed_ptr, 
                                                  score_para_ptr, sn_ratio); 
                                                  */
  if (seed_ptr == nullptr) {
    return false;
  }

  EnvCollPtr env_coll_ptr = env_coll_util::findEnvCollWithSingleEnv(matrix_ptr, seed_ptr, score_para_ptr,
                                                                    sn_ratio); 
  if (env_coll_ptr != nullptr) {
    env_coll_ptr->refineMonoMass(); 
    int score_id = ecscore_list.size();
    ECScorePtr ecscore_ptr = std::make_shared<ECScore>(env_coll_ptr, matrix_ptr, 
                                                       score_id, sn_ratio); 
    ecscore_list.push_back(ecscore_ptr);
    env_coll_ptr->setEcscore(ecscore_ptr->getScore());

    env_coll_list.push_back(env_coll_ptr);
    FracFeaturePtr frac_feature_ptr 
      = env_coll_util::getFracFeature(score_id, ms1_ptr_vec, score_para_ptr->frac_id_, 
                                      score_para_ptr->file_name_, env_coll_ptr, matrix_ptr, 
                                      sn_ratio); 
    frac_feature_ptr->setEcScore(ecscore_ptr->getScore());
    frac_feature_ptr->setHasMs2Spec(true);
    frac_feature_list.push_back(frac_feature_ptr);
    double prec_mono_mass = seed_ptr->getMonoNeutralMass();
    double prec_avg_mass = seed_ptr->getAvgNeutralMass();
    int prec_charge = seed_ptr->getCharge();
    double prec_mono_mz = peak_util::compMz(prec_mono_mass, prec_charge);
    double prec_avg_mz = peak_util::compMz(prec_avg_mass, prec_charge); 
    double prec_inte = seed_ptr->getSeedInte();
    SpecFeaturePtr ms2_feature_ptr = std::make_shared<SpecFeature>(header_ptr, frac_feature_ptr,
                                                                   prec_mono_mz, prec_avg_mz, 
                                                                   prec_charge, prec_inte);
    ms2_feature_list.push_back(ms2_feature_ptr);
    return true;
  }
  else {
    LOG_WARN("Envelope collection is empty!"); 
    return false;
  }
}

void assignEnvColls(FracFeaturePtrVec &frac_feature_list,
                    EnvCollPtrVec &env_coll_list,
                    ECScorePtrVec &ecscore_list,
                    MsMapPtr matrix_ptr,
                    DeconvMsPtrVec &ms1_ptr_vec,
                    MsHeaderPtr2D &ms2_header_ptr_2d,
                    SeedEnvPtr2D &seed_ptr_2d,
                    SpecFeaturePtrVec &ms2_feature_list,
                    TopfdParaPtr topfd_para_ptr,
                    EcscoreParaPtr score_para_ptr) {
  double score_cutoff = topfd_para_ptr->getEcscoreCutoff();
  for (size_t ms1_idx = 0; ms1_idx < ms1_ptr_vec.size(); ms1_idx++) {
    for (size_t i = 0; i < ms2_header_ptr_2d[ms1_idx].size(); i++) {
      MsHeaderPtr header_ptr = ms2_header_ptr_2d[ms1_idx][i];
      bool assigned = getHighestInteEnvColl(frac_feature_list, env_coll_list,  
                                            header_ptr, score_cutoff, ms2_feature_list);
      if (!assigned) {
        // lower the cutoff to 0
        score_cutoff = 0;
        assigned = getHighestInteEnvColl(frac_feature_list, env_coll_list,  
                                         header_ptr, score_cutoff, ms2_feature_list);
      }
      if (!assigned) {
        assigned = getNewEnvColl(header_ptr, ms1_idx, matrix_ptr, score_para_ptr, ecscore_list, 
                                 env_coll_list, ms1_ptr_vec, seed_ptr_2d, 
                                 frac_feature_list, ms2_feature_list); 
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
