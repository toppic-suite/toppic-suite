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

#ifndef TOPPIC_TOPFD_ECSCORE_FEATURE_FEATURE_HPP
#define TOPPIC_TOPFD_ECSCORE_FEATURE_FEATURE_HPP

#include "topfd/ecscore/spectrum/peak_matrix.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"

/*
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/score/get_env_cnn_score.hpp"
#include "topfd/feature_detect/score/get_env_coll_score.hpp"
#include "topfd/feature_detect/score/get_component_score.hpp"
#include "ms/feature/frac_feature.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "topfd/feature_detect/feature_para.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "topfd/feature_detect/env_collection/env_coll_util.hpp"
*/

namespace toppic {

class Feature {
 public:
  Feature(EnvCollPtr env_coll_ptr, PeakMatrixPtr matrix_ptr, int feature_id, double sn_ratio);

  std::vector<float> getEcscoreInput(double max_retention_time);
  
  /*
  Feature(EnvCollPtr env_coll_ptr, PeakMatrixPtr matrix_ptr, int feature_id, double inte);


  Feature(EnvCollection &env_coll, PeakMatrix &peak_matrix, fdeep::model &model, fdeep::model &model_escore,
          int feature_id, double snr);

  FracFeaturePtr static getFeature(int feat_id, DeconvMsPtrVec &ms1_ptr_vec, int frac_id, std::string file_name,
                                   EnvCollection &env_coll, PeakMatrix &peak_matrix, double snr);

  void static
    assign_features(DeconvMsPtrVec &ms1_ptr_vec, const std::string &ms2_file_name,
                    FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                    std::vector<Feature> &features, SpecFeaturePtrVec &ms2_features,
                    std::vector<double> &precMzs, PeakMatrix &peak_matrix, fdeep::model model,
                    const fdeep::model &model_escore, FeatureParaPtr para_ptr, //EnvParaPtr env_para_ptr,
                    TopfdParaPtr topfd_para_ptr);

  bool static
    get_mass_shifted_feature_map(FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                                 FeatureParaPtr para_ptr, MsHeaderPtr hh, double score_thr, double base_mz,
                                 int isolation_windows_mz, SpecFeaturePtrVec &ms2_features);

  bool static
    get_charge_shifted_feature_map(FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                                   FeatureParaPtr para_ptr, MsHeaderPtr hh, double score_thr, double base_mz,
                                   int isolation_windows_mz, SpecFeaturePtrVec &ms2_features);

  bool static
    get_highest_inte_feature_map(FracFeaturePtrVec &frac_features, std::vector<EnvCollection> &env_coll_list,
                                 FeatureParaPtr para_ptr, MsHeaderPtr hh, double score_thr, double base_mz,
                                 int isolation_windows_mz, SpecFeaturePtrVec &ms2_features);

  bool static get_new_feature_map(DeconvMsPtrVec &ms1_ptr_vec, FracFeaturePtrVec &frac_features,
                                  std::vector<EnvCollection> &env_coll_list, std::vector<Feature> &features,
                                  FeatureParaPtr para_ptr, MsHeaderPtr hh, TopfdParaPtr topfd_para_ptr,
                                  SpecFeaturePtrVec &ms2_features, PeakMatrix &peak_matrix, fdeep::model model,
                                  fdeep::model model_escore);

  void static get_empty_feature_map(DeconvMsPtrVec &ms1_ptr_vec, FracFeaturePtrVec &frac_features,
                                    std::vector<EnvCollection> &env_coll_list, std::vector<Feature> &features,
                                    const FeatureParaPtr &para_ptr, MsHeaderPtr hh, TopfdParaPtr topfd_para_ptr,
                                    SpecFeaturePtrVec &ms2_features, PeakMatrix &peak_matrix, fdeep::model model,
                                    fdeep::model model_escore);

  DeconvMsPtrVec static readData(const std::string &file_name);

  double static isMatch(double prec_mass, double feature_mass, const FeatureParaPtr& para_ptr, bool &shift);
  */

  int getFeatureId() const { return feature_id_; }

  void setFeatureId(int feature_id) { feature_id_ = feature_id; }

  int getMinScan() const { return min_scan_; }

  void setMinScan(int min_scan) { min_scan_ = min_scan; }

  int getMaxScan() const { return max_scan_; }

  void setMaxScan(int max_scan) { max_scan_ = max_scan; }

  int getMinCharge() const { return min_charge_; }

  void setMinCharge(int min_charge) { min_charge_ = min_charge; }

  int getMaxCharge() const { return max_charge_; }

  void setMaxCharge(int max_charge) { max_charge_ = max_charge; }

  double getMonoMass() const { return mono_mass_; }

  void setMonoMass(double mono_mass) { mono_mass_ = mono_mass; }

  int getRepCharge() const { return rep_charge_; }

  void setRepCharge(int rep_charge) { rep_charge_ = rep_charge; }

  double getRepMz() const { return rep_mz_; }

  void setRepMz(double rep_mz) { rep_mz_ = rep_mz; }

  double getAbundance() const { return abundance_; }

  void setAbundance(double abundance) { abundance_ = abundance; }

  double getMinElutionTime() const { return min_elution_time_; }

  void setMinElutionTime(double min_elution_time) { min_elution_time_ = min_elution_time; }

  double getMaxElutionTime() const { return max_elution_time_; }

  void setMaxElutionTime(double max_elution_time) { max_elution_time_ = max_elution_time; }

  double getApexElutionTime() const { return apex_elution_time_; }

  void setApexElutionTime(double apex_elution_time) { apex_elution_time_ = apex_elution_time; }

  double getElutionLength() const { return elution_length_; }

  void setElutionLength(double elution_length) { elution_length_ = elution_length; }

  double getMapMaxElutionTime() const { return map_max_elution_time_; }

  void setMapMaxElutionTime(double map_max_elution_time) { map_max_elution_time_ = map_max_elution_time; }

  double getEnvcnnScore() const { return envcnn_score_; }

  void setEnvcnnScore(double envcnn_score) { envcnn_score_ = envcnn_score; }

  double getPercentMatchedPeaks() const { return percent_matched_peaks_; }

  void setPercentMatchedPeaks(double percent_matched_peaks) { percent_matched_peaks_ = percent_matched_peaks; }

  double getIntensityCorrelation() const { return intensity_correlation_; }

  void setIntensityCorrelation(double intensity_correlation) { intensity_correlation_ = intensity_correlation; }

  double getTop3Correlation() const { return top3_correlation_; }

  void setTop3Correlation(double top_3_correlation) { top3_correlation_ = top_3_correlation; }

  double getEvenOddPeakRatio() const { return even_odd_peak_ratio_; }

  void setEvenOddPeakRatio(double even_odd_peak_ratio) { even_odd_peak_ratio_ = even_odd_peak_ratio; }

  double getPercentConsecPeaks() const { return percent_consec_peaks_; }

  void setPercentConsecPeaks(double percent_consec_peaks) { percent_consec_peaks_ = percent_consec_peaks; }

  int getNumTheoPeaks() const { return num_theo_peaks_; }

  void setNumTheoPeaks(int num_theo_peaks) { num_theo_peaks_ = num_theo_peaks; }

  double getMzErrorSum() const { return mz_error_sum_; }

  void setMzErrorSum(double mz_error_sum) { mz_error_sum_ = mz_error_sum; }

  double getScore() const { return score_; }

  void setScore(double score) { score_ = score; }

  int getLabel() const { return label_; }

  void setLabel(int label) { label_ = label; }

  private:
  int feature_id_ = 0;
  int min_scan_ = 0;
  int max_scan_ = 0;
  int min_charge_ = 0;
  int max_charge_ = 0;
  double mono_mass_ = 0;
  int rep_charge_ = 0;
  double rep_mz_ = 0;
  double abundance_ = 0;
  double min_elution_time_ = 0;
  double max_elution_time_ = 0;
  double apex_elution_time_ = 0;
  double elution_length_ = 0;
  double map_max_elution_time_ = 0;
  double envcnn_score_ = 0;
  double percent_matched_peaks_ = 0;
  double intensity_correlation_ = 0;
  double top3_correlation_ = 0;
  double even_odd_peak_ratio_ = 0;
  double percent_consec_peaks_ = 0;
  int num_theo_peaks_ = 0;
  double mz_error_sum_ = 0;
  double score_ = 0;
  int label_ = 0;
};

typedef std::shared_ptr<Feature> FeaturePtr;
typedef std::vector<FeaturePtr> FeaturePtrVec; 

}

#endif 
