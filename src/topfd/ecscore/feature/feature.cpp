//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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
#include "topfd/ecscore/score/component_score.hpp"
#include "topfd/ecscore/score/comp_env_cnn_score.hpp"
#include "topfd/ecscore/score/onnx_ecscore.hpp"
#include "topfd/ecscore/feature/feature.hpp"

namespace toppic {

Feature::Feature(EnvCollPtr env_coll_ptr, MsMapPtr matrix_ptr,
                 int feature_id, double sn_ratio) {
  SeedEnvPtr seed_ptr = env_coll_ptr->getSeedPtr();
  MsMapRowHeaderPtrVec spec_list = matrix_ptr->getHeaderPtrList();

  EnvSetPtr env_set_ptr = env_coll_ptr->getSeedEnvSet();
  feature_id_ = feature_id;

  min_scan_ = env_coll_ptr->getStartSpecId();
  max_scan_ = env_coll_ptr->getEndSpecId();
  min_charge_ = env_coll_ptr->getMinCharge();
  max_charge_ = env_coll_ptr->getMaxCharge();
  mono_mass_ = seed_ptr->getMonoNeutralMass();
  rep_charge_ = seed_ptr->getCharge();
  rep_mz_ = seed_ptr->getMonoMz();

  abundance_ = env_coll_ptr->getIntensity(sn_ratio, matrix_ptr->getBaseInte());

  // convert seconds to minutes
  min_elution_time_ = spec_list[min_scan_]->getRt()/60;
  max_elution_time_ = spec_list[max_scan_]->getRt()/60; 
  int seed_spec_id = seed_ptr->getSpecId();
  apex_elution_time_ = spec_list[seed_spec_id]->getRt()/60;
  elution_length_ = max_elution_time_ - min_elution_time_; 

  double noise_inte = matrix_ptr->getBaseInte();
  EnvSetPtr seed_set_ptr = env_coll_ptr->getSeedEnvSet();
  std::vector<std::vector<double>> theo_map 
    = seed_set_ptr->getScaledTheoIntes(sn_ratio, noise_inte);

  map_max_elution_time_ = spec_list[spec_list.size()-1]->getRt()/60;

  percent_matched_peaks_ = component_score::getMatchedPeakPercent(env_set_ptr, theo_map);
  intensity_correlation_ = component_score::getAggEnvCorr(env_set_ptr);
  top3_correlation_ = component_score::get3ScanCorr(env_set_ptr, seed_spec_id, min_scan_);
  even_odd_peak_ratio_ = component_score::getAggOddEvenPeakRatio(env_set_ptr);
  percent_consec_peaks_ = component_score::getConsecutivePeakPercent(env_set_ptr);
  num_theo_peaks_ = component_score::getTheoPeakNum(theo_map);
  mz_error_sum_ = component_score::getMzErrors(env_set_ptr);

  envcnn_score_ = comp_env_cnn_score::compEnvcnnScore(matrix_ptr, env_coll_ptr); 
  label_ = 0;

  std::vector<float> ecscore_input = getEcscoreInput(map_max_elution_time_);
  score_ = onnx_ecscore::predict(ecscore_input); 
}

std::vector<float> Feature::getEcscoreInput(double max_elution_time) {
  std::vector<float> data;
  data.push_back(envcnn_score_); //1
  data.push_back(elution_length_ / max_elution_time * 2); //2
  data.push_back(percent_matched_peaks_); //3
  data.push_back(rep_charge_); //4
  data.push_back(top3_correlation_); //5
  data.push_back((max_charge_ - min_charge_) / 30.0); //6
  data.push_back(even_odd_peak_ratio_); //7
  return data;
}

void Feature::assignFeatures(FracFeaturePtrVec &frac_feature_list,
                             EnvCollPtrVec &env_coll_list,
                             FeaturePtrVec &feature_list,
                             MsMapPtr matrix_ptr,
                             DeconvMsPtrVec &ms1_ptr_vec,
                             MsHeaderPtr2D &ms2_header_ptr_2d,
                             SeedEnvPtr2D &seed_ptr_2d,
                             SpecFeaturePtrVec &ms2_feature_list,
                             TopfdParaPtr topfd_para_ptr,
                             EcscoreParaPtr score_para_ptr) {
  double score_cutoff = topfd_para_ptr->getEcscoreCutoff();
  for (size_t ms1_id = 0; ms1_id < ms1_ptr_vec.size(); ms1_id++) {
    for (size_t i = 0; i < ms2_header_ptr_2d[ms1_id].size(); i++) {
      MsHeaderPtr header_ptr = ms2_header_ptr_2d[ms1_id][i];
      bool assigned = getHighestInteFeature(frac_feature_list, env_coll_list,  
                                            header_ptr, score_cutoff, ms2_feature_list);
      if (!assigned) {
        // lower the cutoff to 0
        score_cutoff = 0;
        assigned = getHighestInteFeature(frac_feature_list, env_coll_list,  
                                         header_ptr, score_cutoff, ms2_feature_list);
      }
      if (!assigned) {
        assigned = getNewFeature(header_ptr, matrix_ptr, score_para_ptr, feature_list, 
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

bool Feature::getHighestInteFeature(FracFeaturePtrVec &frac_features, EnvCollPtrVec &env_coll_list,
                                    MsHeaderPtr header_ptr, double score_thresh, 
                                    SpecFeaturePtrVec &ms2_features) {
  int ms1_id = header_ptr->getMsOneId();
  double prec_win_bgn = header_ptr->getPrecWinBegin();
  double prec_win_end = header_ptr->getPrecWinEnd();

  SpecFeaturePtrVec new_spec_feats;
  for (size_t coll_id = 0; coll_id < env_coll_list.size(); coll_id++) {
    FracFeaturePtr frac_feature_ptr = frac_features[coll_id];
    // check threshold
    if (frac_feature_ptr->getEcscore() < score_thresh) {
      continue;
    }
    // check retention time range
    if (ms1_id < env_coll_list[coll_id]->getStartSpecId() 
        || ms1_id > env_coll_list[coll_id]->getEndSpecId()) {
      continue;
    }
    EnvSetPtrVec env_sets = env_coll_list[coll_id]->getEnvSetList();
    double feature_mono_mass = env_coll_list[coll_id]->getMass();
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
      std::vector<double> env_intes = env_set_ptr->getXicAllPeakInteList();
      if (env_intes.size() == 0) {
        LOG_WARN("Empty envelope intensity list!");
        continue; 
      }

      double prec_mono_mz = peak_util::compMz(frac_feature_ptr->getMonoMass(), 
                                              env_set_ptr->getCharge());
      int prec_charge = env_set_ptr->getCharge();
      double prec_inte = env_intes[inte_idx];
      if (prec_inte < 0) {
        prec_inte = 0;
      }
      SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(header_ptr,
                                                                 frac_feature_ptr, 
                                                                 prec_mono_mz, 
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

bool Feature::getNewFeature(MsHeaderPtr header_ptr, MsMapPtr matrix_ptr,
                            EcscoreParaPtr score_para_ptr, FeaturePtrVec &feature_list,
                            EnvCollPtrVec &env_coll_list, DeconvMsPtrVec &ms1_ptr_vec,
                            SeedEnvPtr2D &seed_ptr_2d, FracFeaturePtrVec &frac_feature_list,
                            SpecFeaturePtrVec &ms2_feature_list) {
  // set min match envelope to 1 to accept single scan features
  score_para_ptr->min_match_env_ = 1;
  score_para_ptr->min_match_peak_ = 1;
  score_para_ptr->min_seed_match_peak_ = 0;
  double sn_ratio = 0;
  
  int ms1_id = header_ptr->getMsOneId();
  double prec_win_begin = header_ptr->getPrecWinBegin();
  double prec_win_end = header_ptr->getPrecWinEnd();
  SeedEnvPtrVec seed_ptr_list = seed_ptr_2d[ms1_id];
  SeedEnvPtrVec selected_seed_list;
  LOG_DEBUG("ms1 id  " << ms1_id << " seed number " << seed_ptr_list.size());
  for (size_t i = 0; i < seed_ptr_list.size(); i++) {
    double ref_mz = seed_ptr_list[i]->getReferMz();
    if (ref_mz > prec_win_begin && ref_mz < prec_win_end) {
      selected_seed_list.push_back(seed_ptr_list[i]);
    }
  }
  SeedEnvPtr seed_ptr;
  bool valid = false;
  if (selected_seed_list.size() > 0) {
    // choose the highest intensity one
      std::sort(selected_seed_list.begin(), selected_seed_list.end(),
                SeedEnv::cmpSeedInteDec);
    seed_ptr = selected_seed_list[0];  
    valid = seed_env_util::simplePreprocessEnv(matrix_ptr, seed_ptr, 
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
  if (!valid) {
    return false;
  }

  EnvCollPtr env_coll_ptr = env_coll_util::findEnvCollWithSingleEnv(matrix_ptr, seed_ptr, score_para_ptr,
                                                                    sn_ratio); 
  if (env_coll_ptr != nullptr) {
    env_coll_ptr->refineMonoMass(); 
    int feature_id = feature_list.size();
    FeaturePtr feature_ptr = std::make_shared<Feature>(env_coll_ptr, matrix_ptr, 
                                                       feature_id, sn_ratio); 
    feature_list.push_back(feature_ptr);
    env_coll_ptr->setEcscore(feature_ptr->getScore());

    env_coll_list.push_back(env_coll_ptr);
    FracFeaturePtr frac_feature_ptr 
      = env_coll_util::getFracFeature(feature_id, ms1_ptr_vec, score_para_ptr->frac_id_, 
                                      score_para_ptr->file_name_, env_coll_ptr, matrix_ptr, 
                                      sn_ratio); 
    frac_feature_ptr->setEcscore(feature_ptr->getScore());
    frac_feature_ptr->setHasMs2Spec(true);
    frac_feature_list.push_back(frac_feature_ptr);
    double prec_mono_mass = seed_ptr->getMonoNeutralMass();
    int prec_charge = seed_ptr->getCharge();
    double prec_mono_mz = peak_util::compMz(prec_mono_mass, prec_charge);
    double prec_inte = seed_ptr->getSeedInte();
    SpecFeaturePtr ms2_feature_ptr = std::make_shared<SpecFeature>(header_ptr, frac_feature_ptr,
                                                                   prec_mono_mz, prec_charge, prec_inte);
    ms2_feature_list.push_back(ms2_feature_ptr);
    return true;
  }
  else {
    LOG_ERROR("env coll empty"); 
    return false;
  }
}

}
