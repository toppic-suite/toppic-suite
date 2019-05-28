//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_DECONV_FEATURE_FRAC_FEATURE_HPP_
#define TOPPIC_DECONV_FEATURE_FRAC_FEATURE_HPP_

#include <memory>
#include <vector>

namespace toppic {

class FracFeature;
typedef std::shared_ptr<FracFeature> FracFeaturePtr;
typedef std::vector<FracFeaturePtr> FracFeaturePtrVec;

class FracFeature {
 public:
  FracFeature() {}

  FracFeature(int id, int fraction_id, 
              const std::string &file_name,
              double mono_mass, double inte,
              double retent_begin, double retent_end,
              int scan_begin, int scan_end,
              int min_charge, int max_charge);

  FracFeature(std::string line);

  int getId() {return id_;}

  int getFracId() {return frac_id_;}

  std::string getFileName() {return file_name_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return intensity_;}

  double getRetentBegin() {return retent_begin_;}

  double getRetentEnd() {return retent_end_;}

  double getRetentMiddle() {return (retent_begin_ + retent_end_)/2;}

  int getScanBegin() {return scan_begin_;}

  int getScanEnd() {return scan_end_;}

  int getMinCharge() {return min_charge_;}

  int getMaxCharge() {return max_charge_;}

  int getSampleFeatureId() {return sample_feature_id_;}

  double getSampleFeatureInte() {return sample_feature_inte_;}

  void setId(int id) {id_ = id;}

  void setSampleFeatureId(int id) {sample_feature_id_ = id;}

  void setSampleFeatureInte(double inte) {sample_feature_inte_ = inte;}

  static bool cmpMassInc(const FracFeaturePtr &a, const FracFeaturePtr &b) { 
    return a->getMonoMass() < b->getMonoMass();
  }

  static bool cmpInteDec(const FracFeaturePtr &a, const FracFeaturePtr &b) { 
    return a->getIntensity() > b->getIntensity();
  }

  static bool cmpRetentInc(const FracFeaturePtr &a, const FracFeaturePtr &b) { 
    return a->getRetentMiddle() < b->getRetentMiddle();
  }

  static bool cmpFracIncInteDec(const FracFeaturePtr &a, const FracFeaturePtr &b);

 protected:
  int id_;
  int frac_id_;
  std::string file_name_;
  double mono_mass_;
  double intensity_;
  double retent_begin_;
  double retent_end_;
  int scan_begin_;
  int scan_end_;
  int min_charge_;
  int max_charge_;
  int sample_feature_id_ = -1;
  double sample_feature_inte_ = 0;
};

typedef std::vector<FracFeaturePtrVec> FracFeaturePtrVec2D;

}
#endif
