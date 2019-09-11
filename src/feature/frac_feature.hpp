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


#ifndef TOPPIC_FEATURE_FRAC_FEATURE_HPP_
#define TOPPIC_FEATURE_FRAC_FEATURE_HPP_

#include <memory>
#include <vector>

#include "feature/single_charge_feature.hpp"

namespace toppic {

class FracFeature;
typedef std::shared_ptr<FracFeature> FracFeaturePtr;
typedef std::vector<FracFeaturePtr> FracFeaturePtrVec;

class FracFeature {
 public:
  FracFeature() {};

  FracFeature(int id, int fraction_id, 
              const std::string &file_name,
              double mono_mass, double inte,
              int min_ms1_id, int max_ms1_id,
              double retent_begin, double retent_end,
              int scan_begin, int scan_end,
              int min_charge, int max_charge, 
              int env_num);

  FracFeature(std::string line);

  FracFeature(XmlDOMElement* element);

  XmlDOMElement* toXmlElement(XmlDOMDocument* xml_doc);

  int getId() {return id_;}

  int getFracId() {return frac_id_;}

  std::string getFileName() {return file_name_;}

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

  int getEnvNum() {return env_num_;}

  int getSampleFeatureId() {return sample_feature_id_;}

  double getSampleFeatureInte() {return sample_feature_inte_;}

  double getPromexScore() {return promex_score_;}

  bool hasMs2Spec() {return has_ms2_spec_;}

  SingleChargeFeaturePtrVec getSingleFeatures() {return single_features_;}

  void setId(int id) {id_ = id;}

  void setPromexScore(double score) {promex_score_ = score;}

  void setHasMs2Spec(bool has_ms2_spec) {has_ms2_spec_ = has_ms2_spec;}

  void setSampleFeatureId(int id) {sample_feature_id_ = id;}

  void setSampleFeatureInte(double inte) {sample_feature_inte_ = inte;}

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
  int id_;
  int frac_id_;
  std::string file_name_;
  double mono_mass_;
  double intensity_;
  // ms1 ids are used for promex score
  int min_ms1_id_;
  int max_ms1_id_;

  double time_begin_;
  double time_end_;
  int scan_begin_;
  int scan_end_;
  int min_charge_;
  int max_charge_;
  int env_num_ = 0;
  int sample_feature_id_ = -1;
  double sample_feature_inte_ = 0;
  SingleChargeFeaturePtrVec single_features_;

  // used for promex_score
  double promex_score_;
  bool has_ms2_spec_ = false;
};

typedef std::vector<FracFeaturePtrVec> FracFeaturePtrVec2D;

}
#endif
