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

#ifndef TOPPIC_MS_FEATURE_FRAC_FEATURE_HPP_
#define TOPPIC_MS_FEATURE_FRAC_FEATURE_HPP_

#include <memory>
#include <vector>

#include "ms/feature/single_charge_feature.hpp"

namespace toppic {

class FracFeature;
using FracFeaturePtr = std::shared_ptr<FracFeature>;
using FracFeaturePtrVec = std::vector<FracFeaturePtr>;

class FracFeature {
 public:
  FracFeature() {}

  FracFeature(const std::string &file_name,
              int frac_id, int feat_id,  
              double mono_mass, double inte,
              int min_ms1_id, int max_ms1_id,
              double retent_begin, double retent_end,
              int scan_begin, int scan_end,
              int min_charge, int max_charge, 
              double apex_time, int apex_scan, 
              double apex_inte, int rep_charge, 
              double rep_avg_mz, int env_num, 
              double ec_score);  

  explicit FracFeature(const std::string &line);

  explicit FracFeature(XmlDOMElement element);

  XmlDOMElement toXmlElement(XmlDOMDocument* xml_doc, XmlDOMElement parent) const;

  std::string getFileName() const {return file_name_;}

  int getFracId() const {return frac_id_;}

  int getFeatId() const {return feat_id_;}

  double getMonoMass() const {return mono_mass_;}

  double getIntensity() const {return intensity_;}

  int getMinMs1Id() const {return min_ms1_id_;}

  int getMaxMs1Id() const {return max_ms1_id_;}

  double getTimeBegin() const {return time_begin_;}

  double getTimeEnd() const {return time_end_;}

  double getTimeMiddle() const {return (time_begin_ + time_end_)/2;}

  int getScanBegin() const {return scan_begin_;}

  int getScanEnd() const {return scan_end_;}

  int getMinCharge() const {return min_charge_;}

  int getMaxCharge() const {return max_charge_;}

  double getApexTime() const {return apex_time_;}

  int getApexScan() const {return apex_scan_;}

  double getApexInte() const {return apex_inte_;}

  int getRepCharge() const {return rep_charge_;}

  double getRepAvgMz() const {return rep_avg_mz_;}

  int getEnvNum() const {return env_num_;}

  double getEcScore() const {return ec_score_;}

  bool hasMs2Spec() const {return has_ms2_spec_;}

  const SingleChargeFeaturePtrVec& getSingleFeatures() const {return single_features_;}

  void setFracId(int frac_id) {frac_id_ = frac_id;}

  void setFeatId(int feat_id) {feat_id_ = feat_id;}

  void setEcScore(double score) {ec_score_ = score;}

  void setHasMs2Spec(bool has_ms2_spec) {has_ms2_spec_ = has_ms2_spec;}

  void setSingleFeatures(const SingleChargeFeaturePtrVec &single_features) {
    single_features_ = single_features;}

  static bool cmpMassInc(const FracFeaturePtr &a, const FracFeaturePtr &b) { 
    return a->getMonoMass() < b->getMonoMass();
  }

  static bool cmpInteDec(const FracFeaturePtr &a, const FracFeaturePtr &b) { 
    return a->getIntensity() > b->getIntensity();
  }

  static bool cmpTimeInc(const FracFeaturePtr &a, const FracFeaturePtr &b) { 
    return a->getTimeMiddle() < b->getTimeMiddle();
  }

  static bool cmpFracIncInteDec(const FracFeaturePtr &a, const FracFeaturePtr &b);

  static std::string getXmlElementName() {return "frac_feature";}


 protected:
  //mzML file name
  std::string file_name_;
  //the order of the mzML file in the command line input
  int frac_id_;
  //feature id
  int feat_id_;
  double mono_mass_;
  double intensity_;

  // used for ecscore_score
  int min_ms1_id_;
  int max_ms1_id_;

  double time_begin_;
  double time_end_;
  int scan_begin_;
  int scan_end_;
  int min_charge_;
  int max_charge_;
  double apex_time_;
  int apex_scan_;
  double apex_inte_;
  int rep_charge_;
  double rep_avg_mz_;

  int env_num_ = 0;
  double ec_score_;
  bool has_ms2_spec_ = false;

  SingleChargeFeaturePtrVec single_features_;
};

using FracFeaturePtrVec2D = std::vector<FracFeaturePtrVec>;

}  // namespace toppic

#endif
