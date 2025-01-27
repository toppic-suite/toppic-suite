//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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
typedef std::shared_ptr<FracFeature> FracFeaturePtr;
typedef std::vector<FracFeaturePtr> FracFeaturePtrVec;

class FracFeature {
 public:
  FracFeature() {};

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

  FracFeature(std::string line);

  FracFeature(XmlDOMElement* element);

  XmlDOMElement* toXmlElement(XmlDOMDocument* xml_doc);

  std::string getFileName() {return file_name_;}

  int getFracId() {return frac_id_;}

  int getFeatId() {return feat_id_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return intensity_;}

  int getMinMs1Id() {return min_ms1_id_;}

  int getMaxMs1Id() {return max_ms1_id_;}

  double getTimeBegin() {return time_begin_;}

  double getTimeEnd() {return time_end_;}

  double getTimeMiddle() {return (time_begin_ + time_end_)/2;}

  int getScanBegin() {return scan_begin_;}

  int getScanEnd() {return scan_end_;}

  int getMinCharge() {return min_charge_;}

  int getMaxCharge() {return max_charge_;}

  double getApexTime() {return apex_time_;}

  int getApexScan() {return apex_scan_;}

  double getApexInte() {return apex_inte_;}

  int getRepCharge() {return rep_charge_;}

  double getRepAvgMz() {return rep_avg_mz_;}

  int getEnvNum() {return env_num_;}

  double getEcScore() {return ec_score_;}

  bool hasMs2Spec() {return has_ms2_spec_;}

  SingleChargeFeaturePtrVec getSingleFeatures() {return single_features_;}

  void setFracId(int frac_id) {frac_id_ = frac_id;}

  void setFeatId(int feat_id) {feat_id_ = feat_id;}

  void setEcScore(double score) {ec_score_ = score;}

  void setHasMs2Spec(bool has_ms2_spec) {has_ms2_spec_ = has_ms2_spec;}

  void setSingleFeatures(SingleChargeFeaturePtrVec &single_features) {
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

typedef std::vector<FracFeaturePtrVec> FracFeaturePtrVec2D;

}
#endif
