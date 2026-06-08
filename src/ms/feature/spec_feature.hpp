//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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
#include <string>
#include <vector>

#include "ms/spec/peak_util.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/feature/frac_feature.hpp"

namespace toppic {

class SpecFeature;
using SpecFeaturePtr = std::shared_ptr<SpecFeature>;
using SpecFeaturePtrVec = std::vector<SpecFeaturePtr>;

class SpecFeature {
 public:
  explicit SpecFeature(const std::string &line);

  SpecFeature(const MsHeaderPtr &header, const FracFeaturePtr &feature,
              double prec_mono_mz, double prec_avg_mz, 
              int prec_charge, double prec_inte);

  SpecFeature(double prec_mono_mz, double prec_charge); 

  std::string getFileName() const {return file_name_;}

  int getFracId() const {return frac_id_;}

  int getSpecId() const {return spec_id_;}

  std::string getScans() const {return scans_;}

  int getMsOneId() const {return ms_one_id_;}

  int getMsOneScan() const {return ms_one_scan_;}

  double getPrecMonoMz() const {return prec_mono_mz_;}

  double getPrecAvgMz() const {return prec_avg_mz_;}

  double getPrecCharge() const {return prec_charge_;}

  double getPrecMonoMass() const {return peak_util::compPeakNeutralMass(prec_mono_mz_, 
                                                                  prec_charge_);}

  double getPrecInte() const {return prec_inte_;}

  int getFracFeatureId() const {return frac_feature_id_;}

  double getFracFeatureInte() const {return frac_feature_inte_;}

  double getFracFeatureScore() const {return frac_feature_score_;}

  double getFracFeatureMinTime() const {return frac_feature_min_time_;}

  double getFracFeatureMaxTime() const {return frac_feature_max_time_;}

  double getFracFeatureApexTime() const {return frac_feature_apex_time_;}

  void setFracId(int frac_id) {frac_id_ = frac_id;}

  void setSpecId(int id) {spec_id_ = id;}

  void setMsOneId(int id) {ms_one_id_ = id;}

  void setFracFeatureId(int id) {frac_feature_id_ = id;}

  static bool cmpSpecIdInc(const SpecFeaturePtr &a, const SpecFeaturePtr &b) { 
    return a->getSpecId() < b->getSpecId();
  }

  static bool cmpPrecInteDec(const SpecFeaturePtr &a, const SpecFeaturePtr &b) { 
    return a->getPrecInte() > b->getPrecInte();
  }

 protected:
  std::string file_name_;
  int frac_id_;
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
  double prec_mono_mz_;
  double prec_avg_mz_;
  int prec_charge_;
  double prec_inte_;
};

}  // namespace toppic

#endif

