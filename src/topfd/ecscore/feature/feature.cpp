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

#include "common/util/logger.hpp"

#include "ms/spec/peak_util.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/env/env_base.hpp"

#include "topfd/ecscore/score/component_score.hpp"
#include "topfd/ecscore/score/comp_env_cnn_score.hpp"
#include "topfd/ecscore/score/onnx_ecscore.hpp"
#include "topfd/ecscore/envelope/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/spectrum/matrix_spectrum.hpp"
#include "topfd/ecscore/feature/feature.hpp"

namespace toppic {

Feature::Feature(EnvCollPtr env_coll_ptr, PeakMatrixPtr matrix_ptr, 
                 int feature_id, double sn_ratio) {
  SeedEnvelopePtr seed_ptr = env_coll_ptr->getSeedPtr();
  MatrixSpectrumPtrVec spec_list = matrix_ptr->getSpecList();

  EnvSetPtr env_set_ptr = env_coll_ptr->getSeedEnvSet();
  feature_id_ = feature_id;

  min_scan_ = env_coll_ptr->getStartSpecId();
  max_scan_ = env_coll_ptr->getEndSpecId();
  min_charge_ = env_coll_ptr->getMinCharge();
  max_charge_ = env_coll_ptr->getMaxCharge();
  mono_mass_ = seed_ptr->getMass();
  rep_charge_ = seed_ptr->getCharge();
  rep_mz_ = seed_ptr->getPos();

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

void Feature::assignFeatures(const std::string &ms2_file_name, FracFeaturePtrVec &frac_feature_list, 
                             EnvCollPtrVec &env_coll_list, SpecFeaturePtrVec &ms2_feature_list, 
                             std::vector<double> prec_mzs, 
                             TopfdParaPtr topfd_para_ptr, EcscoreParaPtr score_para_ptr, 
                             PeakMatrixPtr matrix_ptr, FeaturePtrVec &feature_list, 
                             DeconvMsPtrVec &ms1_ptr_vec) { 

  double isolation_window_mz = topfd_para_ptr->getPrecWindow();
  double score_cutoff = topfd_para_ptr->getEcscoreCutoff();
  DeconvMsPtrVec ms_ptr_vec = Feature::readData(ms2_file_name);
  //std::cout << "Mapping proteoform features to MS/MS scans started" << std::endl; 
  for (size_t spec_id = 0; spec_id < ms_ptr_vec.size(); spec_id++) {
    MsHeaderPtr header_ptr = ms_ptr_vec[spec_id]->getMsHeaderPtr();
    if (header_ptr->getPrecCharge() == 0) continue;
    double base_mz = prec_mzs[spec_id];
    bool assigned = getHighestInteFeature(frac_feature_list, env_coll_list, score_para_ptr, 
                                          header_ptr, score_cutoff, base_mz, 
                                          isolation_window_mz, ms2_feature_list);
    if (!assigned) {
      // lower the cutoff to 0
      score_cutoff = 0;
      assigned = getHighestInteFeature(frac_feature_list, env_coll_list, score_para_ptr, 
                                       header_ptr, score_cutoff, base_mz, 
                                       isolation_window_mz, ms2_feature_list);
    }
    if (!assigned) {
      assigned = getNewFeature(header_ptr, matrix_ptr, score_para_ptr, feature_list, 
                               env_coll_list, ms1_ptr_vec, frac_feature_list, ms2_feature_list); 
    }
    if (!assigned) {
      LOG_WARN("Scan " << header_ptr->getFirstScanNum() << " does not have MS1 feature!");
    }
  }
  std::sort(ms2_feature_list.begin(), ms2_feature_list.end(), SpecFeature::cmpSpecIdInc);
  MsAlignWriterPtr ms2_ptr = std::make_shared<MsAlignWriter>(ms2_file_name);
  for (const auto &ms_ptr: ms_ptr_vec) {
    ms2_ptr->write(ms_ptr);
  }
  //std::cout << "Mapping proteoform features to MS/MS scans finished" << std::endl; 
}

DeconvMsPtrVec Feature::readData(const std::string &file_name) {
  SimpleMsAlignReader sp_reader(file_name);
  DeconvMsPtrVec ms_ptr_vec;
  DeconvMsPtr ms_ptr;
  while ((ms_ptr = sp_reader.getNextMsPtr()) != nullptr) {
    ms_ptr_vec.push_back(ms_ptr);
  }
  return ms_ptr_vec;
}

bool Feature::getHighestInteFeature(FracFeaturePtrVec &frac_features, EnvCollPtrVec &env_coll_list,
                                    EcscoreParaPtr para_ptr, MsHeaderPtr header_ptr, 
                                    double score_thr, double base_mz,
                                    double isolation_window_mz, SpecFeaturePtrVec &ms2_features) {
  bool assigned = false;
  double env_inte = -1;
  size_t selected_index = -1;
  size_t selected_sub_index = -1;
  int ms1_id = header_ptr->getMsOneId();
  for (size_t feat_id = 0; feat_id < env_coll_list.size(); feat_id++) {
    if (frac_features[feat_id]->getEcscore() < score_thr) {
      continue;
    }
    double feature_mono_mass = env_coll_list[feat_id]->getMass();
    double feature_avg_mass = EnvBase::convertMonoMassToAvgMass(feature_mono_mass); 
    if (ms1_id >= env_coll_list[feat_id]->getStartSpecId() 
        && ms1_id <= env_coll_list[feat_id]->getEndSpecId()) {
      EnvSetPtrVec env_sets = env_coll_list[feat_id]->getEnvSetList();
      for (size_t env_set_id = 0; env_set_id < env_sets.size(); env_set_id++) {
        EnvSetPtr env_set_ptr = env_sets[env_set_id];
        double mz = peak_util::compMz(feature_avg_mass, env_set_ptr->getCharge());
        // this part needs to be improved to consider only peaks in the window
        if ((mz >= base_mz - (isolation_window_mz / 2)) && mz < (base_mz + (isolation_window_mz / 2))) {
          if (ms1_id >= env_set_ptr->getStartSpecId() && env_set_ptr->getEndSpecId()) {
            int inte_idx = ms1_id - env_set_ptr->getStartSpecId();
            std::vector<double> env_intes = env_set_ptr->getXicAllPeakInteList();
            if (env_intes.size() == 0) {
              LOG_WARN("Empty envelope intensity list!");
              continue; 
            }
            double cur_env_inte = env_intes[inte_idx];
            if (env_inte < cur_env_inte) {
              env_inte = cur_env_inte;
              selected_index = feat_id;
              selected_sub_index = env_set_id;
            }
          }
        }
      }
    }
  }
  if (env_inte > -1) {
    FracFeaturePtr feature_ptr = frac_features[selected_index];
    EnvSetPtrVec env_sets = env_coll_list[selected_index]->getEnvSetList();
    EnvSetPtr env_set_ptr = env_sets[selected_sub_index];
    int prec_id = 0;
    double prec_mono_mz = peak_util::compMz(feature_ptr->getMonoMass(), env_set_ptr->getCharge());
    int prec_charge = env_set_ptr->getCharge();
    double prec_inte = env_inte;
    if (prec_inte < 0) {
      prec_inte = 0;
    }
    double apex_time = feature_ptr->getApexTime();
    PrecursorPtr prec_ptr = std::make_shared<Precursor>(prec_id, prec_mono_mz,
                                                        prec_charge, prec_inte,
                                                        apex_time);
    header_ptr->setSinglePrecPtr(prec_ptr);
    SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(header_ptr, feature_ptr);
    ms2_features.push_back(ms2_feature);
    feature_ptr->setHasMs2Spec(true);
    assigned = true;
  }
  return assigned;
}

bool Feature::getNewFeature(MsHeaderPtr header_ptr, PeakMatrixPtr matrix_ptr, 
                            EcscoreParaPtr score_para_ptr, FeaturePtrVec &feature_list, 
                            EnvCollPtrVec &env_coll_list, DeconvMsPtrVec &ms1_ptr_vec, 
                            FracFeaturePtrVec &frac_feature_list, SpecFeaturePtrVec &ms2_feature_list) {
  bool assigned = false;
  // set min match envelope to 1 to accept single scan features
  score_para_ptr->min_match_env_ = 1;
  score_para_ptr->min_match_peak_ = 1;
  double sn_ratio = 0;

  SeedEnvelopePtr seed_ptr = std::make_shared<SeedEnvelope>(header_ptr); 

  bool valid = seed_env_util::simplePreprocessEnv(matrix_ptr, seed_ptr, 
                                                  score_para_ptr, sn_ratio); 
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

    //env_coll.remove_peak_data(peak_matrix); removing peaks is not necessaryhere

    env_coll_list.push_back(env_coll_ptr);
    FracFeaturePtr frac_feature_ptr 
      = env_coll_util::getFracFeature(feature_id, ms1_ptr_vec, score_para_ptr->frac_id_, 
                                      score_para_ptr->file_name_, env_coll_ptr, matrix_ptr, 
                                      sn_ratio); 
    frac_feature_ptr->setEcscore(feature_ptr->getScore());
    frac_feature_ptr->setHasMs2Spec(true);
    frac_feature_list.push_back(frac_feature_ptr);
    SpecFeaturePtr ms2_feature_ptr = std::make_shared<SpecFeature>(header_ptr, frac_feature_ptr);
    ms2_feature_list.push_back(ms2_feature_ptr);
    assigned = true;
  }
  else {
    LOG_ERROR("env coll empty"); 
  }
  return assigned;
}

}
