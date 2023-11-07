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


#ifndef TOPPIC_MS_FEATURE_SPEC_FEATURE_HPP_
#define TOPPIC_MS_FEATURE_SPEC_FEATURE_HPP_

#include <memory>
#include <vector>

#include "ms/spec/peak_util.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/feature/frac_feature.hpp"

namespace toppic {

class SpecFeature;
typedef std::shared_ptr<SpecFeature> SpecFeaturePtr;
typedef std::vector<SpecFeaturePtr> SpecFeaturePtrVec;

class SpecFeature {
 public:
  SpecFeature(std::string line);

  SpecFeature(MsHeaderPtr header, FracFeaturePtr feature,
              double prec_mono_mz, double prec_avg_mz, 
              int prec_charge, double prec_inte);

  SpecFeature(double prec_mono_mz, double prec_charge); 

  int getSpecId() {return spec_id_;}

  int getFracId() {return frac_id_;}

  std::string getFileName() {return file_name_;}
  
  std::string getScans() {return scans_;}

  int getMsOneId() {return ms_one_id_;}

  int getMsOneScan() {return ms_one_scan_;}

  double getPrecMonoMz() {return prec_mono_mz_;}

  double getPrecAvgMz() {return prec_avg_mz_;}

  double getPrecCharge() {return prec_charge_;}

  double getPrecMass() {return peak_util::compPeakNeutralMass(prec_mono_mz_, 
                                                              prec_charge_);}

  double getPrecInte() {return prec_inte_;}

  int getFracFeatureId() {return frac_feature_id_;}

  double getFracFeatureInte() {return frac_feature_inte_;}

  double getFracFeatureScore() {return frac_feature_score_;}

  double getFracFeatureMinTime() {return frac_feature_min_time_;}

  double getFracFeatureMaxTime() {return frac_feature_max_time_;}

  double getFracFeatureApexTime() {return frac_feature_apex_time_;}

  int getSampleFeatureId() {return sample_feature_id_;}

  double getSampleFeatureInte() {return sample_feature_inte_;}

  void setSpecId(int id) {spec_id_ = id;}

  void setMsOneId(int id) {ms_one_id_ = id;}

  void setFracFeatureId(int id) {frac_feature_id_ = id;}

  void setSampleFeatureId(int id) {sample_feature_id_ = id;}

  void setSampleFeatureInte(double inte) {sample_feature_inte_ = inte;}

  static bool cmpSpecIdInc(const SpecFeaturePtr &a, const SpecFeaturePtr &b) { 
    return a->getSpecId() < b->getSpecId();
  }

  static bool cmpPrecInteDec(const SpecFeaturePtr &a, const SpecFeaturePtr &b) { 
    return a->getPrecInte() > b->getPrecInte();
  }

 protected:
  int frac_id_;
  std::string file_name_;
  int spec_id_;
  std::string scans_;
  int ms_one_id_;
  int ms_one_scan_;
  int frac_feature_id_;
  double frac_feature_inte_;
  double frac_feature_score_;
  double frac_feature_min_time_;
  double frac_feature_max_time_;
  double frac_feature_apex_time_;
  int sample_feature_id_;
  double sample_feature_inte_;
  double prec_mono_mz_;
  double prec_avg_mz_;
  int prec_charge_;
  double prec_inte_;
};

}
#endif

