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


#ifndef PROT_FEATURE_FEATURE_HPP_
#define PROT_FEATURE_FEATURE_HPP_

#include <memory>
#include <vector>

namespace prot {

class Feature {
 public:
  Feature(int id, double mono_mass, double inte,
          int scan_begin, int scan_end): 
      id_(id),
      mono_mass_(mono_mass),
      intensity_(inte),
      scan_begin_(scan_begin),
      scan_end_(scan_end) {
      }

  int getId() {return id_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return intensity_;}

  int getScanBegin() {return scan_begin_;}

  int getScanEnd() {return scan_end_;}


 private:
  int id_;
  double mono_mass_;
  double intensity_;
  int scan_begin_;
  int scan_end_;
};

typedef std::shared_ptr<Feature> FeaturePtr;
typedef std::vector<FeaturePtr> FeaturePtrVec;

}
#endif
