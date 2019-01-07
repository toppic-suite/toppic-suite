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


#ifndef TOPPIC_DECONV_FEATURE_FEATURE_HPP_
#define TOPPIC_DECONV_FEATURE_FEATURE_HPP_

#include <memory>
#include <vector>

namespace toppic {

class Feature;
typedef std::shared_ptr<Feature> FeaturePtr;

class Feature {
 public:
  Feature() {}

  Feature(int id, double mono_mass, double inte,
          int scan_begin, int scan_end);

  Feature(int id, double mono_mass, double inte,
          double retent_begin, double retent_end,
          int scan_begin, int scan_end,
          int min_charge, int max_charge);

  int getId() {return id_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return intensity_;}

  double getRetentBegin() {return retent_begin_;}

  double getRetentEnd() {return retent_end_;}

  double getRetentMiddle() {return (retent_begin_ + retent_end_)/2;}

  double getAlignRetentBegin() {return align_retent_begin_;}

  double getAlignRetentEnd() {return align_retent_end_;}

  double getAlignRetentMiddle() {return (align_retent_begin_ + align_retent_end_)/2;}

  int getScanBegin() {return scan_begin_;}

  int getScanEnd() {return scan_end_;}

  int getMinCharge() {return min_charge_;}

  int getMaxCharge() {return max_charge_;}

  int getSampleId() {return sample_id_;}

  void setSampleId(int sample_id) {sample_id_ = sample_id;}

  void setAlignRetentBegin(double begin) {align_retent_begin_ = begin;}

  void setAlignRetentEnd(double end) {align_retent_end_ = end;}

  static bool cmpMassInc(const FeaturePtr &a, const FeaturePtr &b) { 
    return a->getMonoMass() < b->getMonoMass();
  }

  static bool cmpInteDec(const FeaturePtr &a, const FeaturePtr &b) { 
    return a->getIntensity() > b->getIntensity();
  }

  static bool cmpRetentInc(const FeaturePtr &a, const FeaturePtr &b) { 
    return a->getRetentMiddle() < b->getRetentMiddle();
  }

 protected:
  int sample_id_;
  int id_;
  double mono_mass_;
  double intensity_;
  double retent_begin_;
  double retent_end_;
  double align_retent_begin_;
  double align_retent_end_;
  int scan_begin_;
  int scan_end_;
  int min_charge_;
  int max_charge_;
};

typedef std::vector<FeaturePtr> FeaturePtrVec;

}
#endif
