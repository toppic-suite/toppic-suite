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


#ifndef TOPPIC_FEATURE_SAMPLE_FEATURE_HPP_
#define TOPPIC_FEATURE_SAMPLE_FEATURE_HPP_

#include <memory>
#include <vector>

#include "feature/frac_feature.hpp"

namespace toppic {

class SampleFeature;
typedef std::shared_ptr<SampleFeature> SampleFeaturePtr;

class SampleFeature {
 public:
  SampleFeature() {}

  int getSampleId() {return sample_id_;}

  int getId() {return id_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return intensity_;}

  double getRetentBegin() {return retent_begin_;}

  double getRetentEnd() {return retent_end_;}

  double getRetentMiddle() {return (retent_begin_ + retent_end_)/2;}

  int getMinCharge() {return min_charge_;}

  int getMaxCharge() {return max_charge_;}

  void setSampleId(int sample_id) {sample_id_ = sample_id;}

  static bool cmpMassInc(const SampleFeaturePtr &a, const SampleFeaturePtr &b) { 
    return a->getMonoMass() < b->getMonoMass();
  }

  static bool cmpInteDec(const SampleFeaturePtr &a, const SampleFeaturePtr &b) { 
    return a->getIntensity() > b->getIntensity();
  }

  static bool cmpRetentInc(const SampleFeaturePtr &a, const SampleFeaturePtr &b) { 
    return a->getRetentMiddle() < b->getRetentMiddle();
  }

 protected:
  int sample_id_;
  int id_;
  double mono_mass_;
  double intensity_;
  double retent_begin_;
  double retent_end_;
  int min_charge_;
  int max_charge_;
  FracFeaturePtrVec frac_features_;
};

typedef std::vector<SampleFeaturePtr> SampleFeaturePtrVec;

}
#endif
