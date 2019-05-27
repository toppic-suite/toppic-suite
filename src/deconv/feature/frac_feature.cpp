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

#include <fstream>

#include "common/util/str_util.hpp"
#include "deconv/feature/frac_feature.hpp"

namespace toppic {

FracFeature::FracFeature(int id, int fraction_id, double mono_mass, double inte,
                         double retent_begin, double retent_end,
                         int scan_begin, int scan_end,
                         int min_charge, int max_charge): 
    id_(id),
    fraction_id_(fraction_id),
    mono_mass_(mono_mass),
    intensity_(inte),
    retent_begin_(retent_begin),
    retent_end_(retent_end),
    scan_begin_(scan_begin),
    scan_end_(scan_end),
    min_charge_(min_charge),
    max_charge_(max_charge) {
    }

FracFeature::FracFeature(std::string line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  id_ = std::stoi(strs[0]);
  fraction_id_ = std::stoi(strs[1]);
  mono_mass_ = std::stod(strs[2]);
  intensity_ = std::stod(strs[3]);
  retent_begin_ = std::stod(strs[4]);
  retent_end_ = std::stod(strs[5]);
  scan_begin_ = std::stoi(strs[6]);
  scan_end_ = std::stoi(strs[7]);
  min_charge_ = std::stoi(strs[8]);
  max_charge_ = std::stoi(strs[9]);
}

void FracFeature::writeFeatures(const std::string &output_file_name,
                                const FracFeaturePtrVec &features) {
  std::ofstream of(output_file_name);
  of.precision(16);
  of << "ID" << "\t"
      << "Fraction ID" << "\t"
      << "Mass" << "\t"
      << "Intensity" << "\t"
      << "Time begin" << "\t"
      << "Time end" << "\t"
      << "First scan" << "\t"
      << "Last scan" << "\t"
      << "Minimum charge state" << "\t"
      << "Maximum charge state" 
      << std::endl;
  for (size_t i = 0; i < features.size(); i++) {
    FracFeaturePtr feature = features[i];
    of << feature->getId() << "\t"
        << feature->getFractionId() << "\t"
        << feature->getMonoMass() << "\t"
        << feature->getIntensity() << "\t"
        << feature->getRetentBegin() << "\t"
        << feature->getRetentEnd() << "\t"
        << feature->getScanBegin() << "\t"
        << feature->getScanEnd() << "\t"
        << feature->getMinCharge() << "\t"
        << feature->getMaxCharge() << "\t"
        << std::endl;
  }
  of.close();
}

}
