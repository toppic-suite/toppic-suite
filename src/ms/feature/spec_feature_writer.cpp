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

#include <set>
#include <iomanip> 
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/feature/spec_feature_writer.hpp"

namespace toppic {

namespace spec_feature_writer {

void writeHeader(std::ofstream &of) {
  of << "File_name" << "\t"
     << "Fraction_ID" << "\t"
     << "Spectrum_ID" << "\t"
     << "Scans" << "\t"
     << "MS_one_ID" << "\t"
     << "MS_one_scans" << "\t"
     << "Fraction_feature_ID" << "\t"
     << "Fraction_feature_intensity" << "\t"
     << "Fraction_feature_score" << "\t"
     << "Fraction_feature_min_time" << "\t"
     << "Fraction_feature_max_time" << "\t"
     << "Fraction_feature_apex_time" << "\t"
     << "Precursor_monoisotopic_mz" << "\t"
     << "Precursor_average_mz" << "\t"
     << "Precursor_charge" << "\t"
     << "Precursor_intensity"
     << std::endl;
}

void writeOneFeature(std::ofstream &of, SpecFeaturePtr feature) {
  of << feature->getFileName() << "\t"
     << feature->getFracId() << "\t"
     << feature->getSpecId() << "\t"
     << feature->getScans() << "\t"
     << feature->getMsOneId() << "\t"
     << feature->getMsOneScan() << "\t"
     << feature->getFracFeatureId() << "\t"
     << feature->getFracFeatureInte() << "\t"
     << feature->getFracFeatureScore() << "\t"
     << feature->getFracFeatureMinTime()/60 << "\t"
     << feature->getFracFeatureMaxTime()/60 << "\t"
     << feature->getFracFeatureApexTime()/60 << "\t"
     << std::fixed << std::setprecision(5) 
     << feature->getPrecMonoMz() << "\t"
     << feature->getPrecAvgMz() << "\t"
     << std::setprecision(0) 
     << feature->getPrecCharge() << "\t"
     << std::setprecision(5) 
     << feature->getPrecInte() 
     << std::endl;
}

void writeFeatures(const std::string &output_file_name,
                   const SpecFeaturePtrVec &features) {
  std::ofstream of(output_file_name);
  writeHeader(of);

  for (size_t i = 0; i < features.size(); i++) {
    SpecFeaturePtr feature = features[i];
    writeOneFeature(of, feature);
  }
  of.close();
}

}

} /* namespace toppic */
