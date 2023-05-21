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

#include "topfd/ecscore/score/component_score.hpp"
#include "topfd/ecscore/score/comp_env_cnn_score.hpp"
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

  double scan_max_rt = spec_list[spec_list.size()-1]->getRt();
  min_elution_time_ = spec_list[min_scan_]->getRt()/scan_max_rt;
  max_elution_time_ = spec_list[max_scan_]->getRt()/scan_max_rt; 
  int seed_spec_id = seed_ptr->getSpecId();
  apex_elution_time_ = spec_list[seed_spec_id]->getRt()/scan_max_rt;
  elution_length_ = max_elution_time_ - min_elution_time_; 

  double noise_inte = matrix_ptr->getBaseInte();
  EnvSetPtr seed_set_ptr = env_coll_ptr->getSeedEnvSet();
  std::vector<std::vector<double>> theo_map 
    = seed_set_ptr->getScaledTheoIntes(sn_ratio, noise_inte);

  percent_matched_peaks_ = component_score::getMatchedPeakPercent(env_set_ptr, theo_map);
  intensity_correlation_ = component_score::getAggEnvCorr(env_set_ptr);
  top3_correlation_ = component_score::get3ScanCorr(env_set_ptr, seed_spec_id, min_scan_);
  even_odd_peak_ratio_ = component_score::getAggOddEvenPeakRatio(env_set_ptr);
  percent_consec_peaks_ = component_score::getConsecutivePeakPercent(env_set_ptr);
  num_theo_peaks_ = component_score::getTheoPeakNum(theo_map);
  mz_error_sum_ = component_score::getMzErrors(env_set_ptr);

  envcnn_score_ = comp_env_cnn_score::compEnvcnnScore(matrix_ptr, env_coll_ptr); 
  label_ = 0;

  /*
  std::vector<double> data;
  data.push_back(envcnn_score_); //1
  data.push_back(elution_length_ / 60.0); //2
  data.push_back(percent_matched_peaks_); //3
  data.push_back(std::log(abundance_)); //4
  data.push_back(rep_charge_); //5
  data.push_back(top3_correlation_); //6
  data.push_back((max_charge_ - min_charge_) / 30.0); //7
  data.push_back(even_odd_peak_ratios_); //8
  score_ = env_coll_score::get_env_coll_score(model_escore, data);
  */
}

/*
Feature::Feature(EnvCollPtr env_coll_ptr, PeakMatrixPtr matrix_ptr, 
                 int feature_id, double inte) {
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
  abundance_ = inte;

  double scan_max_rt = spec_list[spec_list.size()-1]->getRt();

  min_elution_time_ = spec_list[min_scan_]->getRt() / scan_max_rt;
  max_elution_time_ = spec_list[max_scan_]->getRt() / scan_max_rt; 
  int seed_scan = seed_ptr->getSpecId();
  apex_elution_time_ = spec_list[seed_scan]->getRt() / scan_max_rt;
  elution_length_ = (max_elution_time_ - min_elution_time_) /scan_max_rt; 

  percent_matched_peaks_ = 0;
  intensity_correlation_ = 0;
  top3_correlation_ = 0;
  even_odd_peak_ratio_ = 0;
  percent_consec_peaks_ = 0;
  num_theo_peaks_ = 0;
  mz_error_sum_ = 0;
  envcnn_score_ = 0;
  score_ = 0;
  label_ = 0;
}
*/

/*
DeconvMsPtrVec Feature::readData(const std::string &file_name) {
  SimpleMsAlignReader sp_reader(file_name);
  DeconvMsPtrVec ms_ptr_vec;
  DeconvMsPtr ms_ptr;
  while ((ms_ptr = sp_reader.getNextMsPtr()) != nullptr) {
    ms_ptr_vec.push_back(ms_ptr);
  }
  return ms_ptr_vec;
}

FracFeaturePtr Feature::getFeature(int feat_id, DeconvMsPtrVec &ms1_ptr_vec, int frac_id, std::string file_name,
                                   EnvCollection &env_coll, PeakMatrix &peak_matrix, double snr) {
  double noise_inte = peak_matrix.get_min_inte();
  spec_list spectra_list = peak_matrix.get_spectra_list();
  int ms1_id_begin = env_coll.getStartSpecId();
  int ms1_id_end = env_coll.getEndSpecId();
  double feat_inte = env_coll.get_intensity(snr, peak_matrix.get_min_inte());
  double feat_mass = env_coll.getMass();
  int min_charge = env_coll.getMinCharge();
  int max_charge = env_coll.getMaxCharge();
  double ms1_time_begin = env_coll.get_min_elution_time(spectra_list);
  double ms1_time_end = env_coll.get_max_elution_time(spectra_list);
  int ms1_scan_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getFirstScanNum();
  int ms1_scan_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getFirstScanNum();

  // get apex inte
  double time_apex = env_coll.get_apex_elution_time(spectra_list);
  EnvSet es = env_coll.get_seed_env_set();
  std::vector<double> env_intes = es.getEnvIntes();
  double apex_inte = 0;
  if (env_intes.size() > 0) {
    int inte_idx = env_coll.getBaseSpecID() - env_coll.getStartSpecId();
    apex_inte = env_intes[inte_idx];
  }

  FracFeaturePtr feature_ptr = std::make_shared<FracFeature>(feat_id, frac_id, file_name, feat_mass, feat_inte,
                                                             ms1_id_begin, ms1_id_end, ms1_time_begin, ms1_time_end,
                                                             ms1_scan_begin, ms1_scan_end, min_charge, max_charge,
                                                             0, time_apex, apex_inte);
  SingleChargeFeaturePtrVec single_features;
  for (EnvSet &es: env_coll.getEnvSetList()) {
    int id_begin = es.getStartSpecId();
    int id_end = es.getEndSpecId();
    double time_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getRetentionTime();
    double time_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getRetentionTime();
    int scan_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getFirstScanNum();
    int scan_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getFirstScanNum();
    double inte = es.comp_intensity(snr, noise_inte);
    int charge = es.getCharge();
    std::vector<double> xic = es.getXicEnvIntes();
    SingleChargeFeaturePtr single_feature = std::make_shared<SingleChargeFeature>(charge, time_begin, time_end,
                                                                                  scan_begin, scan_end,
                                                                                  inte, 0, id_begin, id_end,
                                                                                  feat_mass, xic);
    single_features.push_back(single_feature);
  }
  feature_ptr->setSingleFeatures(single_features);
  return feature_ptr;
}

double Feature::isMatch(double prec_mass, double feature_mass, const FeatureParaPtr &para_ptr, bool &shift) {
  std::vector<double> search_masses = para_ptr->getSearchMasses(prec_mass);
  double min_diff = std::numeric_limits<double>::max();
  for (size_t j = 0; j < search_masses.size(); j++) {
    double mass_diff = std::abs(search_masses[j] - feature_mass);
    if (mass_diff < min_diff) {
      if (j > 0)
        shift = true;
      min_diff = mass_diff;
    }
  }
  return min_diff;
}

void Feature::assign_features(DeconvMsPtrVec &ms1_ptr_vec, const std::string &ms2_file_name,
                              FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                              std::vector<Feature> &features, SpecFeaturePtrVec &ms2_features,
                              std::vector<double> &precMzs, PeakMatrix &peak_matrix, fdeep::model model,
                              const fdeep::model &model_escore, FeatureParaPtr para_ptr, TopfdParaPtr topfd_para_ptr) {

  int isolation_windows_mz = topfd_para_ptr->getPrecWindow();
  DeconvMsPtrVec ms_ptr_vec = Feature::readData(ms2_file_name);
  std::cout << "\r" << "Mapping Proteoforms Features on MS2 Scans" << std::flush;
  for (size_t spec_id = 0; spec_id < ms_ptr_vec.size(); spec_id++) {
    MsHeaderPtr hh = ms_ptr_vec[spec_id]->getMsHeaderPtr();
    if (hh->getPrecCharge() == 0) continue;

    double base_mz = precMzs[spec_id];
    bool assigned_status = get_mass_shifted_feature_map(frac_features, env_coll_list, para_ptr, hh, topfd_para_ptr->getECScore(),
                                                        base_mz, isolation_windows_mz, ms2_features);
    if (assigned_status) continue;

    assigned_status = get_charge_shifted_feature_map(frac_features, env_coll_list, para_ptr, hh, topfd_para_ptr->getECScore(), base_mz,
                                                     isolation_windows_mz, ms2_features);
    if (assigned_status) continue;

    assigned_status = get_highest_inte_feature_map(frac_features, env_coll_list, para_ptr, hh, topfd_para_ptr->getECScore(), base_mz,
                                                   isolation_windows_mz, ms2_features);
    if (assigned_status) continue;

    assigned_status = get_new_feature_map(ms1_ptr_vec, frac_features, env_coll_list, features, para_ptr, hh, topfd_para_ptr,
                                          ms2_features, peak_matrix, model, model_escore);
    if (assigned_status) continue;

    get_empty_feature_map(ms1_ptr_vec, frac_features, env_coll_list, features, para_ptr, hh, topfd_para_ptr,
                          ms2_features, peak_matrix, model, model_escore);
    if (assigned_status) continue;

  }
  std::sort(ms2_features.begin(), ms2_features.end(), SpecFeature::cmpSpecIdInc);
  MsAlignWriterPtr ms2_ptr = std::make_shared<MsAlignWriter>(ms2_file_name);
  for (const auto &ms_ptr: ms_ptr_vec)
    ms2_ptr->write(ms_ptr);
}

void Feature::get_empty_feature_map(DeconvMsPtrVec &ms1_ptr_vec, FracFeaturePtrVec &frac_features,
                                    std::vector<EnvCollection> &env_coll_list, std::vector<Feature> &features,
                                    const FeatureParaPtr &para_ptr, MsHeaderPtr hh, TopfdParaPtr topfd_para_ptr,
                                    SpecFeaturePtrVec &ms2_features, PeakMatrix &peak_matrix, fdeep::model model,
                                    fdeep::model model_escore) {
  int env_coll_num = env_coll_list.size();
  SeedEnvelope env = SeedEnvelope(hh);
  std::vector<ExpEnvelope> env_list;
  std::vector<EnvSet> env_set_list;
  int min_charge = env.getCharge();
  int max_charge = env.getCharge();
  int start_spec_id = env.getSpecId();
  int end_spec_id = env.getSpecId();
  EnvSet es = EnvSet(env, env_list, start_spec_id, end_spec_id, peak_matrix.get_min_inte(),
                     topfd_para_ptr->getMsOneSnRatio());
  env_set_list.push_back(es);

  EnvCollection env_coll = EnvCollection(env, env_set_list, min_charge, max_charge, start_spec_id, end_spec_id);
  Feature feature = Feature(env_coll, peak_matrix, env_coll_num, env.getInte());

  features.push_back(feature);
  env_coll.setEcscore(feature.getScore());
  env_coll.remove_peak_data(peak_matrix);
  env_coll_list.push_back(env_coll);
  FracFeaturePtr feature_ptr = getFeature(env_coll_num, ms1_ptr_vec, para_ptr->frac_id_, para_ptr->file_name_,
                                          env_coll, peak_matrix, topfd_para_ptr->getMsOneSnRatio());
  feature_ptr->setPromexScore(feature.getScore());
  frac_features.push_back(feature_ptr);
  SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(hh, feature_ptr);
  ms2_features.push_back(ms2_feature);
  feature_ptr->setHasMs2Spec(true);
}

bool Feature::get_new_feature_map(DeconvMsPtrVec &ms1_ptr_vec, FracFeaturePtrVec &frac_features,
                                  std::vector<EnvCollection> &env_coll_list, std::vector<Feature> &features,
                                  FeatureParaPtr para_ptr, MsHeaderPtr hh, TopfdParaPtr topfd_para_ptr,
                                  SpecFeaturePtrVec &ms2_features, PeakMatrix &peak_matrix, fdeep::model model,
                                  fdeep::model model_escore) {
  bool assigned_status = false;
  para_ptr->match_peak_tole_ = 0; // Set to 0 to include single scan features

  SeedEnvelope env = SeedEnvelope(hh);
  double min_mz = peak_matrix.get_min_mz() - para_ptr->mass_tole_;
  double max_mz = peak_matrix.get_max_mz() + para_ptr->mass_tole_;
  env.rm_peaks(min_mz, max_mz);
  env_set_util::comp_peak_start_end_idx(peak_matrix, env, para_ptr->mass_tole_);
  EnvCollection env_coll = env_coll_util::find_env_collection(peak_matrix, env, para_ptr,
                                                              topfd_para_ptr->getMsOneSnRatio());
  if (!env_coll.isEmpty()) {
    int env_coll_num = env_coll_list.size();
    env_coll.refine_mono_mass();
    Feature feature = Feature(env_coll, peak_matrix, model, model_escore, env_coll_num,
                              topfd_para_ptr->getMsOneSnRatio());
    features.push_back(feature);
    env_coll.setEcscore(feature.getScore());
    env_coll.remove_peak_data(peak_matrix);
    env_coll_list.push_back(env_coll);
    FracFeaturePtr feature_ptr = getFeature(env_coll_num, ms1_ptr_vec, para_ptr->frac_id_, para_ptr->file_name_,
                                            env_coll, peak_matrix, topfd_para_ptr->getMsOneSnRatio());
    feature_ptr->setPromexScore(feature.getScore());
    frac_features.push_back(feature_ptr);
    SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(hh, feature_ptr);
    ms2_features.push_back(ms2_feature);
    feature_ptr->setHasMs2Spec(true);
    assigned_status = true;
  }
  return assigned_status;
}

bool
Feature::get_highest_inte_feature_map(FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                                      FeatureParaPtr para_ptr, MsHeaderPtr hh, double score_thr, double base_mz,
                                      int isolation_windows_mz, SpecFeaturePtrVec &ms2_features) {
  bool assigned_status = false;
  double env_inte = -1;
  size_t selected_index = -1;
  size_t selected_sub_index = -1;
  int ms1_id = hh->getMsOneId();
  for (size_t feat_id = 0; feat_id < env_coll_list.size(); feat_id++) {
    if (frac_features[feat_id]->getPromexScore() < score_thr)
      continue;
    double feature_mass = env_coll_list[feat_id].getMass();
    if (ms1_id >= env_coll_list[feat_id].getStartSpecId() and ms1_id <= env_coll_list[feat_id].getEndSpecId()) {
      std::vector<EnvSet> sfs = env_coll_list[feat_id].getEnvSetList();
      for (size_t sf_id = 0; sf_id < sfs.size(); sf_id++) {
        EnvSet sf = sfs[sf_id];
        double mz = get_mz(feature_mass, sf.getCharge());
        if ((mz >= base_mz - (isolation_windows_mz / 2)) and mz < (base_mz + (isolation_windows_mz / 2))) {
          if (ms1_id >= sf.getStartSpecId() and ms1_id <= sf.getEndSpecId()) {
            int inte_idx = ms1_id - sf.getStartSpecId();
            std::vector<double> env_intes = sf.getEnvIntes();
            if (env_intes.size() ==
                0) //////////////////////////////////////////////////////////////////////////////////////// ERRORRRR
              continue; //////////////////////////////////////////////////////////////////////////////////////// ERRORRRR
            double sf_env_inte = env_intes[inte_idx];
            if (env_inte < sf_env_inte) {
              env_inte = sf_env_inte;
              selected_index = feat_id;
              selected_sub_index = sf_id;
            }
          }
        }
      }
    }
  }
  if (env_inte > -1) {
    FracFeaturePtr feature_ptr = frac_features[selected_index];
    std::vector<EnvSet> sfs = env_coll_list[selected_index].getEnvSetList();
    EnvSet shortlisted_single_Charge_feature = sfs[selected_sub_index];
    hh->setPrecMonoMz(get_mz(feature_ptr->getMonoMass(), shortlisted_single_Charge_feature.getCharge()));
    hh->setPrecCharge(shortlisted_single_Charge_feature.getCharge());
    if (env_inte > 0)
      hh->setPrecInte(env_inte);
    SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(hh, feature_ptr);
    ms2_features.push_back(ms2_feature);
    feature_ptr->setHasMs2Spec(true);
    assigned_status = true;
  }
  return assigned_status;
}

bool
Feature::get_charge_shifted_feature_map(FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                                        FeatureParaPtr para_ptr, MsHeaderPtr hh, double score_thr, double base_mz,
                                        int isolation_windows_mz, SpecFeaturePtrVec &ms2_features) {
  bool assigned_status = false;
  double env_inte = -1;
  size_t selected_index = -1;
  size_t selected_sub_index = -1;
  int ms1_id = hh->getMsOneId();
  for (size_t feat_id = 0; feat_id < env_coll_list.size(); feat_id++) {
    if (frac_features[feat_id]->getPromexScore() < score_thr)
      continue;
    bool shift = false;
    double feature_mass = env_coll_list[feat_id].getMass();
    double prec_mz = hh->getPrecMonoMz();
    double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mz);
    if (ms1_id >= env_coll_list[feat_id].getStartSpecId() and ms1_id <= env_coll_list[feat_id].getEndSpecId()) {
      std::vector<EnvSet> sfs = env_coll_list[feat_id].getEnvSetList();
      for (size_t sf_id = 0; sf_id < sfs.size(); sf_id++) {
        EnvSet sf = sfs[sf_id];
        double mz = get_mz(feature_mass, sf.getCharge());
        if ((mz >= base_mz - (isolation_windows_mz / 2)) and mz < (base_mz + (isolation_windows_mz / 2))) {
          if (ms1_id >= sf.getStartSpecId() and ms1_id <= sf.getEndSpecId()) {
            int inte_idx = ms1_id - sf.getStartSpecId();
            std::vector<double> env_intes = sf.getEnvIntes();
            if (env_intes.size() ==
                0) //////////////////////////////////////////////////////////////////////////////////////// ERRORRRR
              continue; //////////////////////////////////////////////////////////////////////////////////////// ERRORRRR
            double sf_env_inte = env_intes[inte_idx];
            double diff = isMatch(prec_mz, mz, para_ptr, shift);
            if (diff > error_tole)
              continue;
            if (env_inte < sf_env_inte) {
              env_inte = sf_env_inte;
              selected_index = feat_id;
              selected_sub_index = sf_id;
            }
          }
        }
      }
    }
  }
  if (env_inte > -1) {
    FracFeaturePtr feature_ptr = frac_features[selected_index];
    std::vector<EnvSet> sfs = env_coll_list[selected_index].getEnvSetList();
    EnvSet shortlisted_single_Charge_feature = sfs[selected_sub_index];
    hh->setPrecMonoMz(get_mz(feature_ptr->getMonoMass(), shortlisted_single_Charge_feature.getCharge()));
    hh->setPrecCharge(shortlisted_single_Charge_feature.getCharge());
    if (env_inte > 0)
      hh->setPrecInte(env_inte);
    SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(hh, feature_ptr);
    ms2_features.push_back(ms2_feature);
    feature_ptr->setHasMs2Spec(true);
    assigned_status = true;
  }
  return assigned_status;
}

bool
Feature::get_mass_shifted_feature_map(FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                                      FeatureParaPtr para_ptr, MsHeaderPtr hh, double score_thr, double base_mz,
                                      int isolation_windows_mz, SpecFeaturePtrVec &ms2_features) {
  bool assigned_status = false;
  double mass_diff = std::numeric_limits<double>::max();
  double env_inte = -1;
  size_t selected_index = -1;
  size_t selected_sub_index = -1;
  bool selected_shift = false;
  int ms1_id = hh->getMsOneId();
  for (size_t feat_id = 0; feat_id < env_coll_list.size(); feat_id++) {
    if (frac_features[feat_id]->getPromexScore() < score_thr)
      continue;
    bool shift = false;
    double feature_mass = env_coll_list[feat_id].getMass();
    double prec_mass = hh->getPrecMonoMass();
    double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
    double diff = isMatch(prec_mass, feature_mass, para_ptr, shift);
    if (diff > error_tole)
      continue;
    if (ms1_id >= env_coll_list[feat_id].getStartSpecId() and ms1_id <= env_coll_list[feat_id].getEndSpecId()) {
      std::vector<EnvSet> sfs = env_coll_list[feat_id].getEnvSetList();
      for (size_t sf_id = 0; sf_id < sfs.size(); sf_id++) {
        EnvSet sf = sfs[sf_id];
        double mz = get_mz(feature_mass, sf.getCharge());
        if ((mz >= base_mz - (isolation_windows_mz / 2)) and mz < (base_mz + (isolation_windows_mz / 2))) {
          if (ms1_id >= sf.getStartSpecId() and ms1_id <= sf.getEndSpecId()) {
            int inte_idx = ms1_id - sf.getStartSpecId();
            std::vector<double> env_intes = sf.getEnvIntes();
            if (env_intes.size() ==
                0) //////////////////////////////////////////////////////////////////////////////////////// ERRORRRR
              continue; //////////////////////////////////////////////////////////////////////////////////////// ERRORRRR
            double sf_env_inte = env_intes[inte_idx];
            if (mass_diff > diff) {
              mass_diff = diff;
              env_inte = sf_env_inte;
              selected_index = feat_id;
              selected_sub_index = sf_id;
              selected_shift = shift;
            }
          }
        }
      }
    }
  }
  if (env_inte > -1) {
    FracFeaturePtr feature_ptr = frac_features[selected_index];
    std::vector<EnvSet> sfs = env_coll_list[selected_index].getEnvSetList();
    EnvSet shortlisted_single_Charge_feature = sfs[selected_sub_index];
    if (selected_shift) {
      hh->setPrecMonoMz(get_mz(feature_ptr->getMonoMass(), shortlisted_single_Charge_feature.getCharge()));
      hh->setPrecCharge(shortlisted_single_Charge_feature.getCharge());
    }
    hh->setPrecInte(env_inte);
    SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(hh, feature_ptr);
    ms2_features.push_back(ms2_feature);
    feature_ptr->setHasMs2Spec(true);
    assigned_status = true;
  }
  return assigned_status;
}
*/
}
