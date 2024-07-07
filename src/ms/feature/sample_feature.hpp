//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_MS_FEATURE_SAMPLE_FEATURE_HPP_
#define TOPPIC_MS_FEATURE_SAMPLE_FEATURE_HPP_

#include <memory>
#include <string>
#include <vector>

#include "ms/feature/frac_feature.hpp"

namespace toppic {

class SampleFeature;
typedef std::shared_ptr<SampleFeature> SampleFeaturePtr;

class SampleFeature {
 public:
  SampleFeature() {}

  SampleFeature(const std::string &line);

  SampleFeature(FracFeaturePtr frac_feature, int id);

  SampleFeature(FracFeaturePtrVec &frac_features, int id);

  void init(FracFeaturePtr frac_feature); 

  int getSampleId() {return sample_id_;}

  int getFeatId() {return feat_id_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return intensity_;}

  double getTimeBegin() {return time_begin_;}

  double getTimeEnd() {return time_end_;}

  double getElutionLen() {return time_end_ - time_begin_;}

  double getMinScan() {return min_scan_;}

  double getMaxScan() {return max_scan_;}

  int getMinCharge() {return min_charge_;}

  int getMaxCharge() {return max_charge_;}

  double getTimeMiddle() {return (time_begin_ + time_end_)/2;}

  int getRepFracId() {return rep_frac_id_;}

  double getRepApexTime() {return rep_apex_time_;}

  double getRepApexScan() {return rep_apex_scan_;}

  double getRepCharge() {return rep_charge_;}

  double getRepAvgMz() {return rep_avg_mz_;}

  double getRepApexInte() {return rep_apex_inte_;}

  int getEnvNum() {return env_num_;}

  double getRepEcScore() {return rep_ec_score_;}

  int getMinFracId() {return min_frac_id_;}

  int getMaxFracId() {return max_frac_id_;}

  void setSampleId(int sample_id) {sample_id_ = sample_id;}

  static bool cmpMassInc(const SampleFeaturePtr &a, const SampleFeaturePtr &b) { 
    return a->getMonoMass() < b->getMonoMass();
  }

  static bool cmpInteDec(const SampleFeaturePtr &a, const SampleFeaturePtr &b) { 
    return a->getIntensity() > b->getIntensity();
  }

  static bool cmpTimeInc(const SampleFeaturePtr &a, const SampleFeaturePtr &b) { 
    return a->getTimeMiddle() < b->getTimeMiddle();
  }

 protected:
  int sample_id_ = 0;
  int feat_id_;
  double mono_mass_;
  double intensity_;
  double time_begin_;
  double time_end_;
  int min_charge_;
  int max_charge_;
  int min_frac_id_;
  int max_frac_id_;
  int min_scan_;
  int max_scan_;

  int rep_frac_id_;
  double rep_apex_time_;
  int rep_apex_scan_;
  int rep_charge_;
  double rep_avg_mz_;
  double rep_apex_inte_;
  int env_num_;
  double rep_ec_score_;
  //FracFeaturePtrVec frac_features_;
};

typedef std::vector<SampleFeaturePtr> SampleFeaturePtrVec;

}
#endif
