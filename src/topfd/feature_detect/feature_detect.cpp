//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <cmath>
#include <iostream>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/time_util.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/mod_util.hpp"
#include "seq/fasta_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/match_env.hpp"
#include "ms/feature/frac_feature.hpp"
#include "ms/feature/single_charge_feature.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/feature/peak_cluster.hpp"
#include "ms/feature/sample_feature.hpp"
#include "ms/feature/sample_feature_writer.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "topfd/feature_detect/feature_para.hpp"
#include "topfd/feature_detect/feature_detect.hpp"

namespace toppic {

namespace feature_detect {

bool peakExists (DeconvMsPtrVec &ms1_ptr_vec, DeconvPeakPtr peak) {
  int sp_id = peak->getSpId();
  int peak_id = peak->getId();
  DeconvPeakPtr remain_peak = ms1_ptr_vec[sp_id]->getPeakPtr(peak_id);
  if (remain_peak == nullptr) {
    return false;
  }
  else {
    return true;
  }
}

void getMatchedPeaks(DeconvMsPtrVec &ms1_ptr_vec, double prec_mass,
                     DeconvPeakPtrVec &matched_peaks, 
                     int ms1_id_begin, int ms1_id_end, 
                     FeatureParaPtr para_ptr) {
  std::vector<double> mass_diff_vec;
  std::vector<double> error_tole_vec;  
  if (ms1_ptr_vec.size() == 0) return;
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  std::vector<double> ext_masses = para_ptr->getExtendMasses(prec_mass);
  for (int i = ms1_id_begin; i <= ms1_id_end; i++) {
    DeconvMsPtr ms1_ptr = ms1_ptr_vec[i];
    for (size_t j = 0; j < ms1_ptr->size(); j++) {
      DeconvPeakPtr peak = ms1_ptr->getPeakPtr(j);
      if (peak != nullptr) {
        for (size_t k = 0; k < ext_masses.size(); k++) {
          double mass_diff = std::abs(ext_masses[k] - peak->getPosition());
          if (mass_diff <= error_tole) {
            peak->setSpId(i);
            matched_peaks.push_back(peak);
            break;
          }
          else{
            mass_diff_vec.push_back(mass_diff);
            error_tole_vec.push_back(error_tole);
          }
        }
      }
    }
  }
  if (matched_peaks.size() > 0){}
  else{
    for (int i = 0; i < 100; i++){
      std::cout << "mass diff: " << mass_diff_vec[i] << ", error_tole : " << error_tole_vec[i] << std::endl;
    }
  }
}

bool containPrecursor(DeconvMsPtr ms1_ptr, double prec_mass, int charge, FeatureParaPtr para_ptr) {
  if (ms1_ptr == nullptr) return false;
  // double prec_chrg = best_ptr->getPrecCharge();
  std::vector<double> ext_masses = para_ptr->getExtendMasses(prec_mass);

  double min_diff = std::numeric_limits<double>::max();
  for (size_t i = 0; i < ms1_ptr->size(); i++) {
    DeconvPeakPtr peak = ms1_ptr->getPeakPtr(i);
    // if (peak->getCharge() == prec_chrg) {
    // do not test charge
    if (peak != nullptr) {
      if (charge < 0 || peak->getCharge() == charge) {
        for (size_t j = 0; j < ext_masses.size(); j++) {
          double mass_diff = std::abs(ext_masses[j] - peak->getPosition());
          if (mass_diff < min_diff) {
            min_diff = mass_diff;
          }
        }
      }
    }
  }

  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  if (min_diff <= error_tole) {
    return true;
  }
  return false;
}

int getMs1IdBegin(const DeconvMsPtrVec &ms1_ptr_vec, int sp_id, double prec_mass,
                  FeatureParaPtr para_ptr) {
  int cur_id = sp_id;
  int result_id = sp_id;
  int miss_num = 0;
  while (cur_id > 0) {
    if (containPrecursor(ms1_ptr_vec[cur_id], prec_mass, -1, para_ptr)) {
      miss_num = 0;
      result_id = cur_id;
    } else {
      miss_num++;
    }
    cur_id--;
    if (miss_num >= 15) {
      break;
    }
  }
  return result_id;
}

int getMs1IdEnd(const DeconvMsPtrVec &ms1_ptr_vec, int sp_id, double prec_mass,
                FeatureParaPtr para_ptr) {
  int cur_id = sp_id;
  int result_id = sp_id;
  int miss_num = 0;
  while (cur_id < static_cast<int>(ms1_ptr_vec.size()) ) {
    if (containPrecursor(ms1_ptr_vec[cur_id], prec_mass, -1, para_ptr)) {
      miss_num = 0;
      result_id = cur_id;
    } else {
      miss_num++;
    }
    cur_id++;
    if (miss_num >= 15) {
      break;
    }
  }
  return result_id;
}

std::pair<int,int> getSingleMs1Range(const DeconvMsPtrVec &ms1_ptr_vec, double prec_mass, 
                                     int charge, int begin_id, int end_id, FeatureParaPtr para_ptr) {
  int begin = end_id;
  int end = begin_id;
  for (int i = begin_id; i <= end_id; i++) {
    if (containPrecursor(ms1_ptr_vec[i], prec_mass, charge, para_ptr)) {
      if (i < begin) {
        begin = i;
      }
      if (i > end) {
        end = i;
      }
    }
  }
  std::pair<int,int> result(begin, end); 
  return result;
}


double getFeatureInte(DeconvPeakPtrVec &matched_peaks) {
  double inte = 0;
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    inte += matched_peaks[i]->getIntensity();
  }
  return inte;
}

double getFeatureMass(double best_peak_mass, DeconvPeakPtrVec &matched_peaks,
                      FeatureParaPtr para_ptr) {
  std::vector<double> offsets = para_ptr->getExtendOffsets();
  size_t offset_num = offsets.size();
  DeconvPeakPtrVec2D offset_peaks;
  for (size_t i = 0; i < offsets.size(); i++) {
    DeconvPeakPtrVec peaks;
    offset_peaks.push_back(peaks);
  }
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(best_peak_mass);
  // get peaks for each offset
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    DeconvPeakPtr peak = matched_peaks[i];
    for (size_t j = 0; j < offset_num; j++) {
      double mass_diff = std::abs(best_peak_mass + offsets[j] - peak->getPosition());
      if (mass_diff <= error_tole) {
        offset_peaks[j].push_back(peak);
        break;
      }
    }
  }
  // get feature intensities for each offset
  std::vector<double> intes(offset_num, 0);
  for (size_t i = 0; i < offset_num; i++) {
    for (size_t j = 0; j < offset_peaks[i].size(); j++) {
      intes[i] += offset_peaks[i][j]->getIntensity();
    }
  }
  // get the best offset
  int best_offset = -1;
  double best_inte = -1;
  for (size_t i = 0; i < offset_num; i++) {
    if (intes[i] > best_inte) {
      best_inte = intes[i];
      best_offset = i;
    }
  }
  //get the mass as weighted average
  double inte_sum = 0;
  double product_sum = 0;
  for (size_t i = 0; i < offset_peaks[best_offset].size(); i++) {
    DeconvPeakPtr peak = offset_peaks[best_offset][i];
    inte_sum += peak->getIntensity();
    product_sum = product_sum + peak->getPosition() * peak->getIntensity();
  }
  if (inte_sum == 0.0) {
    LOG_ERROR("Intensity sum is 0!");
  }
  return product_sum/inte_sum;
}

int getMinCharge (DeconvPeakPtrVec &matched_peaks) {
  int min_charge = std::numeric_limits<int>::max();
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    if (matched_peaks[i]->getCharge() < min_charge) {
      min_charge = matched_peaks[i]->getCharge();
    }
  }
  return min_charge;
}

int getMaxCharge (DeconvPeakPtrVec &matched_peaks) {
  int max_charge = -1;
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    if (matched_peaks[i]->getCharge() > max_charge) {
      max_charge = matched_peaks[i]->getCharge();
    }
  }
  return max_charge;
}

DeconvPeakPtrVec getPeaksWithCharge(DeconvPeakPtrVec &matched_peaks, int charge) {
  DeconvPeakPtrVec peaks; 
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    if (matched_peaks[i]->getCharge() == charge) {
      peaks.push_back(matched_peaks[i]);
    }
  }
  return peaks;
}

FracFeaturePtr getFeature(int sp_id, double prec_mass, int feat_id, 
                          DeconvMsPtrVec &ms1_ptr_vec,
                          DeconvPeakPtrVec &matched_peaks, FeatureParaPtr para_ptr) {
  int ms1_id_begin = getMs1IdBegin(ms1_ptr_vec, sp_id, prec_mass, para_ptr);
  int ms1_id_end = getMs1IdEnd(ms1_ptr_vec, sp_id, prec_mass, para_ptr);
  getMatchedPeaks(ms1_ptr_vec, prec_mass, matched_peaks,
                  ms1_id_begin, ms1_id_end, para_ptr);
  if (matched_peaks.size() == 0) {
    return nullptr;
  }
  double feat_inte = getFeatureInte(matched_peaks);
  double feat_mass = getFeatureMass(prec_mass, matched_peaks, para_ptr);
  int min_charge = getMinCharge(matched_peaks);
  int max_charge = getMaxCharge(matched_peaks);
  double ms1_time_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getRetentionTime();
  double ms1_time_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getRetentionTime();
  int ms1_scan_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getFirstScanNum();
  int ms1_scan_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getFirstScanNum();
  FracFeaturePtr feature_ptr = std::make_shared<FracFeature>(feat_id, 
                                                             para_ptr->frac_id_,
                                                             para_ptr->file_name_,
                                                             feat_mass,
                                                             feat_inte,
                                                             ms1_id_begin,
                                                             ms1_id_end,
                                                             ms1_time_begin,
                                                             ms1_time_end,
                                                             ms1_scan_begin,
                                                             ms1_scan_end, 
                                                             min_charge,
                                                             max_charge, 
                                                             matched_peaks.size());
  SingleChargeFeaturePtrVec single_features;
  for (int charge = min_charge; charge <= max_charge; charge++) {
    DeconvPeakPtrVec peaks = getPeaksWithCharge(matched_peaks, charge);
    if (peaks.size() > 0) {
      std::pair<int, int> range = getSingleMs1Range(ms1_ptr_vec, prec_mass, charge, ms1_id_begin, 
                                                    ms1_id_end, para_ptr);
      int id_begin = range.first;
      int id_end = range.second;
      //LOG_ERROR("Id begin" << id_begin << " Id end " << id_end);
      double time_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getRetentionTime();
      double time_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getRetentionTime();
      int scan_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getFirstScanNum();
      int scan_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getFirstScanNum();
      double inte = getFeatureInte(peaks);
      int peak_num = peaks.size();
      SingleChargeFeaturePtr single_feature = std::make_shared<SingleChargeFeature>(charge,
                                                                                    time_begin, time_end,
                                                                                    scan_begin, scan_end,
                                                                                    inte, peak_num);
      single_features.push_back(single_feature);
    }
  }
  feature_ptr->setSingleFeatures(single_features);
  return feature_ptr;
}

void removePeaks (DeconvMsPtrVec &ms1_ptr_vec, DeconvPeakPtrVec &matched_peaks) {
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    DeconvPeakPtr peak = matched_peaks[i];
    int sp_id = peak->getSpId();
    int peak_id = peak->getId();
    ms1_ptr_vec[sp_id]->setPeakPtrNull(peak_id);
  }
}

MatchEnvPtr getMatchEnv(const PeakPtrVec &peak_list, int sp_id, double mono_neutral_mass, int charge,
                        EnvParaPtr env_para_ptr) {

  EnvelopePtr ref_env = EnvBase::getStaticEnvByMonoMass(mono_neutral_mass);
  if (ref_env == nullptr) {
    LOG_ERROR("reference envelope is null");
    exit(EXIT_FAILURE);
  }
  double new_mono_charge_mass = mono_neutral_mass + charge * mass_constant::getIsotopeMass();
  // original envelope has charge 1
  double ori_mono_charge_mass = ref_env->getMonoNeutralMass() + mass_constant::getIsotopeMass();
  double mass_diff = new_mono_charge_mass - ori_mono_charge_mass;
  EnvelopePtr theo_env = ref_env->convertToTheo(mass_diff, charge); 

  int mass_group = env_para_ptr->getMassGroup(mono_neutral_mass);
  // LOG_DEBUG("theo env raw complete");

  // get the highest 85%--95% peaks
  double percentage = env_para_ptr->getPercentBound(mass_group);
  // include all possible peaks;
  double min_inte = 0.0;
  int max_back_peak_num = 1000;
  int max_forw_peak_num = 1000;
  theo_env = theo_env->getSubEnv(percentage, min_inte, max_back_peak_num, max_forw_peak_num); 

  RealEnvPtr real_env = std::make_shared<RealEnv>(peak_list, theo_env, env_para_ptr->getMzTolerance(),
                                                  min_inte);
  if (real_env == nullptr) {
    LOG_ERROR("real env is null");
  }
  real_env->setSpId(sp_id);
  //LOG_DEBUG("theo peak_num " << theo_env->getPeakNum() << " real peak num " << real_env->getPeakNum());
  MatchEnvPtr match_env = std::make_shared<MatchEnv>(mass_group, theo_env, real_env);
  return match_env;
}

void findMsOneFeatures(DeconvMsPtrVec &ms1_ptr_vec, PeakPtrVec2D & raw_peaks, 
                       FeatureParaPtr para_ptr, FracFeaturePtrVec &features,
                       EnvParaPtr env_para_ptr) {
  //get all peaks
  DeconvPeakPtrVec all_peaks;
  for (size_t i = 0; i < ms1_ptr_vec.size(); i++) {
    DeconvPeakPtrVec peaks = ms1_ptr_vec[i]->getPeakPtrVec();
    all_peaks.insert(std::end(all_peaks), std::begin(peaks), std::end(peaks));
  }
  //sort all peaks by intensities
  std::sort(all_peaks.begin(), all_peaks.end(), Peak::cmpInteDec);
  int feat_id = 0;
  size_t peak_idx = 0;
  while (feat_id < para_ptr->feature_num_ && peak_idx < all_peaks.size()) {
    DeconvPeakPtr best_peak = all_peaks[peak_idx];
    if (peakExists(ms1_ptr_vec, best_peak)) {
      //std::cout << "Find feature " << feat_id << " peak intensity " << best_peak->getIntensity() << std::endl; 
      int ref_sp_id = best_peak->getSpId();
      double prec_mass = best_peak->getPosition();
      DeconvPeakPtrVec matched_peaks;
      LOG_DEBUG("feature id " << feat_id);
      FracFeaturePtr feature_ptr = getFeature(ref_sp_id, prec_mass, feat_id, ms1_ptr_vec,
                                              matched_peaks, para_ptr);
      if (feature_ptr == nullptr) {
        LOG_ERROR("Empty feature!");
        exit(EXIT_FAILURE);
      }
      else{
        // check if the feature has at least 2 envelopes
        //if (feature_ptr->getEnvNum() > 1) {
          double ref_mono_mass = feature_ptr->getMonoMass();
          double ref_charge = best_peak->getCharge();
          int sp_id = best_peak->getSpId();
          MatchEnvPtr match_env = getMatchEnv(raw_peaks[sp_id], sp_id, ref_mono_mass, ref_charge, env_para_ptr); 
          if (match_env == nullptr) {
            LOG_ERROR("matche env is null");
          }
          PeakClusterPtr peak_cluster = std::make_shared<PeakCluster>(match_env->getTheoEnvPtr());
          RealEnvPtrVec real_envs; 
          for (size_t i = 0; i < matched_peaks.size(); i++) {
            ref_charge = matched_peaks[i]->getCharge();
            sp_id = matched_peaks[i]->getSpId();
            match_env = getMatchEnv(raw_peaks[sp_id], sp_id, ref_mono_mass, ref_charge, env_para_ptr);
            real_envs.push_back(match_env->getRealEnvPtr());
          }
          LOG_DEBUG("get real envs done");
          peak_cluster->addEnvelopes(feature_ptr, real_envs); 
          LOG_DEBUG("add real envs done");
          bool check_pvalue = true;
          peak_cluster->updateScore(raw_peaks, check_pvalue);
          LOG_DEBUG("update score done");
          double promex_score = para_ptr->peak_cluster_score_ptr_->getScore(peak_cluster);
          LOG_DEBUG("get promex score done");
          feature_ptr->setPromexScore(promex_score);
          features.push_back(feature_ptr);
          
        //}
        removePeaks(ms1_ptr_vec, matched_peaks);
      }
      feat_id++;
    }
    peak_idx++;
  }
}

void readHeaders(const std::string & file_name, MsHeaderPtrVec &header_ptr_vec) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(file_name, sp_num_in_group, nullptr,
                          std::set<std::string>());
  DeconvMsPtr ms_ptr;
  LOG_DEBUG("Start search");
  while ((ms_ptr = sp_reader.getNextMs()) != nullptr) {
    header_ptr_vec.push_back(ms_ptr->getMsHeaderPtr());
    //std::cout << std::flush <<  "reading spectrum " << header_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  //std::cout << std::endl;
}

inline bool isMatch(FracFeaturePtr feature_ptr, MsHeaderPtr header, FeatureParaPtr para_ptr) {
  int ms1_scan = header->getMsOneScan();
  if (ms1_scan < feature_ptr->getScanBegin()) {
    return false;
  }
  if (ms1_scan > feature_ptr->getScanEnd()) {
    return false;
  }
  double prec_mass = header->getPrecMonoMass();
  std::vector<double> search_masses = para_ptr->getSearchMasses(prec_mass);

  double feature_mass = feature_ptr->getMonoMass();

  double min_diff = std::numeric_limits<double>::max();
  for (size_t j = 0; j < search_masses.size(); j++) {
    double mass_diff = std::abs(search_masses[j] - feature_mass);
    if (mass_diff < min_diff) {
      min_diff = mass_diff;
    }
  }
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  if (min_diff <= error_tole) {
    return true;
  }
  return false;
}

inline FracFeaturePtr getMatchedFeaturePtr(FracFeaturePtrVec &features, 
                                           MsHeaderPtr header,
                                           FeatureParaPtr para_ptr) {
  for (size_t i = 0; i < features.size(); i++) {
    if (isMatch(features[i], header, para_ptr)) {
      return features[i];
    }
  }
  return nullptr;
}

void getMs2Features(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtrVec &header_ptr_vec, 
                    FracFeaturePtrVec &features, FeatureParaPtr para_ptr, 
                    SpecFeaturePtrVec &ms2_features) {
  MsHeaderPtrVec sorted_ptrs = header_ptr_vec;
  std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), MsHeader::cmpPrecInteDec);
  for (size_t i = 0; i < sorted_ptrs.size(); i++) {
    MsHeaderPtr header = sorted_ptrs[i];
    FracFeaturePtr ft_ptr = getMatchedFeaturePtr(features, header, para_ptr);
    if (ft_ptr != nullptr) {
      SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(header, ft_ptr);
      ms2_features.push_back(ms2_feature);
      ft_ptr->setHasMs2Spec(true);
    }
    else {
      int sp_id = header->getMsOneId();
      double prec_mass = header->getPrecMonoMass();
      if (prec_mass > 0) {
        int feat_id = static_cast<int>(features.size());

        DeconvPeakPtrVec matched_peaks;
        FracFeaturePtr feature_ptr = getFeature(sp_id, prec_mass, feat_id, ms1_ptr_vec,
                                                matched_peaks, para_ptr);
        // if we find a feature in ms1.msalign
        // it is possible that some ms headers do not have matched features. 
        if (feature_ptr != nullptr) {
          feature_ptr->setHasMs2Spec(true);
          feature_ptr->setPromexScore(-1000);
          features.push_back(feature_ptr);
          removePeaks(ms1_ptr_vec, matched_peaks);

          SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(header, feature_ptr);
          ms2_features.push_back(ms2_feature);
        }
        else {
          LOG_WARN("Cannot find features in LC/MS! Spectrum id: " << sp_id);
        }
      }
    }
  }
  std::sort(ms2_features.begin(), ms2_features.end(), SpecFeature::cmpSpecIdInc);
}

void getSampleFeatures(SampleFeaturePtrVec &sample_features, FracFeaturePtrVec &frac_features,
                       SpecFeaturePtrVec &spec_features) {
  //sample features;
  for (size_t i = 0; i < frac_features.size(); i++) {
    SampleFeaturePtr sample_feature = std::make_shared<SampleFeature>(frac_features[i], frac_features[i]->getId());
    sample_features.push_back(sample_feature);
    frac_features[i]->setSampleFeatureId(frac_features[i]->getId());
    frac_features[i]->setSampleFeatureInte(frac_features[i]->getIntensity());
  }

  //spec features
  std::map<int,FracFeaturePtr> feature_map;
  for (size_t i = 0; i < frac_features.size(); i++) {
    feature_map[frac_features[i]->getId()] =  frac_features[i];
  }

  for (size_t j = 0; j < spec_features.size(); j++) {
    SpecFeaturePtr spec_feature = spec_features[j];
    int frac_feature_id = spec_feature->getFracFeatureId(); 
    FracFeaturePtr frac_feature = feature_map.find(frac_feature_id)->second;
    spec_feature->setSampleFeatureId(frac_feature->getSampleFeatureId());
    spec_feature->setSampleFeatureInte(frac_feature->getSampleFeatureInte());
  }
}

void process(int frac_id, const std::string &sp_file_name, 
             bool missing_level_one, const std::string &resource_dir) { 
  //logger::setLogLevel(2);
  FeatureParaPtr para_ptr 
      = std::make_shared<FeaturePara>(frac_id, sp_file_name, resource_dir);
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
  std::string base_name = file_util::basename(sp_file_name);
  // read ms1 deconvoluted spectra
  DeconvMsPtrVec ms1_ptr_vec;
  FracFeaturePtrVec frac_features;
  if (!missing_level_one) {
    std::string ms1_file_name = base_name + "_ms1.msalign";
    MsAlignReader::readMsOneSpectra(ms1_file_name, ms1_ptr_vec);
    PeakPtrVec2D raw_peaks; 
    RawMsReaderPtr raw_reader_ptr = std::make_shared<RawMsReader>(sp_file_name);
    raw_reader_ptr->getMs1Peaks(raw_peaks);
    raw_reader_ptr = nullptr;
    findMsOneFeatures(ms1_ptr_vec, raw_peaks, para_ptr, frac_features, env_para_ptr);
  
  
  LOG_DEBUG("start reading ms2");
  std::string ms2_file_name = base_name + "_ms2.msalign";
  MsHeaderPtrVec header_ptr_vec;
  readHeaders(ms2_file_name, header_ptr_vec);

  SpecFeaturePtrVec ms2_features;
  getMs2Features(ms1_ptr_vec, header_ptr_vec, frac_features, para_ptr, ms2_features);

  SampleFeaturePtrVec sample_features;
  getSampleFeatures(sample_features, frac_features, ms2_features);

  std::string output_file_name = base_name + "_feature.xml";
  frac_feature_writer::writeXmlFeatures(output_file_name, frac_features);
  std::string batmass_file_name = base_name + "_frac.mzrt.csv";
  frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);
  std::string sample_feature_file_name = base_name + "_ms1.feature";
  sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);

  output_file_name = base_name + "_ms2.feature";
  spec_feature_writer::writeFeatures(output_file_name, ms2_features); 
  }
}

}  // namespace 

}  // namespace toppic 

