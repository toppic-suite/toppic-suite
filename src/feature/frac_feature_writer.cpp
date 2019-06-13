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
#include "spec/peak.hpp"
#include "deconv/env/envelope.hpp"
#include "deconv/env/env_base.hpp"
#include "feature/frac_feature_writer.hpp"

namespace toppic {

namespace frac_feature_writer {

void writeHeader(std::ofstream &of) {
  of.precision(16);
  of << "ID" << "\t"
      << "Fraction_ID" << "\t"
      << "File_name" << "\t"
      << "Mass" << "\t"
      << "Intensity" << "\t"
      << "Time_begin" << "\t"
      << "Time_end" << "\t"
      << "First_scan" << "\t"
      << "Last_scan" << "\t"
      << "Minimum_charge_state" << "\t"
      << "Maximum_charge_state" << "\t"
      << "Sample_feature_Id" << "\t"
      << "Sample_feature_intensity"
      << std::endl;
}

void writeOneFeature(std::ofstream &of, FracFeaturePtr feature) {
  of << feature->getId() << "\t"
      << feature->getFracId() << "\t"
      << feature->getFileName() << "\t"
      << feature->getMonoMass() << "\t"
      << feature->getIntensity() << "\t"
      << feature->getTimeBegin() << "\t"
      << feature->getTimeEnd() << "\t"
      << feature->getScanBegin() << "\t"
      << feature->getScanEnd() << "\t"
      << feature->getMinCharge() << "\t"
      << feature->getMaxCharge() << "\t"
      << feature->getSampleFeatureId() << "\t"
      << feature->getSampleFeatureInte() 
      << std::endl;
}

void writeFeatures(const std::string &output_file_name,
                   const FracFeaturePtrVec &features) {
  std::ofstream of(output_file_name);
  writeHeader(of);

  for (size_t i = 0; i < features.size(); i++) {
    FracFeaturePtr feature = features[i];
    writeOneFeature(of, feature);
  }
  of.close();
}

void writeBatMassFeatures(const std::string &output_file_name,
                          const FracFeaturePtrVec &features) {
  std::ofstream of(output_file_name);
  std::string delimit = ",";
  of << "ID" << delimit
      << "Fraction_ID" << delimit
      << "Mass" << delimit
      << "MonoMz" << delimit
      << "Charge" << delimit
      << "Intensity" << delimit
      << "mzLo" << delimit
      << "mzHi" << delimit
      << "rtLo" << delimit
      << "rtHi" << delimit
      << "color" << delimit
      << "opacity" << delimit
      << "Sample_feature_Id" << delimit
      << "Sample_feature_intensity"
      << std::endl;
  for (size_t i = 0; i < features.size(); i++) {
    FracFeaturePtr feature = features[i];
    int min_charge = feature->getMinCharge();
    int max_charge = feature->getMaxCharge();
    for (int j = min_charge; j <= max_charge; j++) {
      double mono_mass = feature->getMonoMass(); 
      double mono_mz = Peak::compMonoMz(mono_mass, j);
      EnvelopePtr ref_env = EnvBase::getStaticEnvByMonoMass(mono_mass);
      EnvelopePtr theo_env = ref_env->distrToTheoMono(mono_mz, j);
      double min_inte = 0.03;
      EnvelopePtr filtered_env = theo_env->getSubEnv(min_inte); 
      of << feature->getId() << delimit
          << feature->getFracId() << delimit
          << feature->getMonoMass() << delimit
          << mono_mz << delimit
          << j << delimit
          << feature->getIntensity() << delimit
          << filtered_env->getMinMz() << delimit
          << filtered_env->getMaxMz() << delimit
          << (feature->getTimeBegin()/60) << delimit
          << (feature->getTimeEnd()/60) << delimit
          << "#FF0000" << delimit
          << "0.1" << delimit
          << feature->getSampleFeatureId() << delimit
          << feature->getSampleFeatureInte() 
          << std::endl;
    }
  }
  of.close();
}

}

} /* namespace toppic */
