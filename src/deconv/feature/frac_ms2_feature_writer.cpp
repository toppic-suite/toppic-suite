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

#include <set>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "deconv/feature/frac_ms2_feature_writer.hpp"

namespace toppic {

namespace frac_ms2_feature_writer {

  int id_;
  int frac_id_;
  std::string file_name_;
  std::string scans_;
  int ms_one_id_;
  std::string ms_one_scans_;
  double prec_mass_;
  double prec_inte_;
  int frac_feature_id_;
  double frac_feature_inte_;
  int sample_feature_id_;
  double sample_feature_inte_;

void writeHeader(std::ofstream &of) {
  of.precision(16);
  of << "ID" << "\t"
      << "Fraction_ID" << "\t"
      << "File_name" << "\t"
      << "Scans" << "\t"
      << "MS_one_ID" << "\t"
      << "MS_one_scans" << "\t"
      << "Precursor_mass" << "\t"
      << "Precursor_intensity" << "\t"
      << "Fraction_feature_ID" << "\t"
      << "Fraction_feature_intensity" << "\t"
      << "Sample_feature_ID" << "\t"
      << "Sample_feature_intensity"
      << std::endl;
}

void writeOneFeature(std::ofstream &of, FracMs2FeaturePtr feature) {
  of << feature->getId() << "\t"
      << feature->getFracId() << "\t"
      << feature->getFileName() << "\t"
      << feature->getScans() << "\t"
      << feature->getMsOneId() << "\t"
      << feature->getMsOneScan() << "\t"
      << feature->getPrecMass() << "\t"
      << feature->getPrecInte() << "\t"
      << feature->getFracFeatureId() << "\t"
      << feature->getFracFeatureInte() << "\t"
      << feature->getSampleFeatureId() << "\t"
      << feature->getSampleFeatureInte() 
      << std::endl;
}

void writeFeatures(const std::string &output_file_name,
                   const FracMs2FeaturePtrVec &features) {
  std::ofstream of(output_file_name);
  writeHeader(of);

  for (size_t i = 0; i < features.size(); i++) {
    FracMs2FeaturePtr feature = features[i];
    writeOneFeature(of, feature);
  }
  of.close();
}

}

} /* namespace toppic */
