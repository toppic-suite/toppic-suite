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

#include <set>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_impl.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "ms/spec/peak.hpp"
#include "ms/env/envelope.hpp"
#include "ms/env/env_base.hpp"
#include "ms/feature/frac_feature_writer.hpp"

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
      << "Envelope_num" << "\t"
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
      << feature->getEnvNum() << "\t"
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
      << "Envelope_num" << delimit
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
      << "promex_score" 
      << std::endl;
  for (size_t i = 0; i < features.size(); i++) {
    FracFeaturePtr feature = features[i];
    double mono_mass = feature->getMonoMass(); 
    SingleChargeFeaturePtrVec single_features = feature->getSingleFeatures();

    for (size_t j = 0; j < single_features.size(); j++) {
      SingleChargeFeaturePtr single_feature = single_features[j];
      int charge = single_feature->getCharge();
      double mono_mz = Peak::compMz(mono_mass, charge);
      EnvelopePtr ref_env = EnvBase::getStaticEnvByMonoMass(mono_mass);
      EnvelopePtr theo_env = ref_env->distrToTheoMono(mono_mz, charge);
      double min_inte = 0.03;
      EnvelopePtr filtered_env = theo_env->getSubEnv(min_inte); 
      //margin for envelopes
      double margin = 0.1; 
      double min_mz = filtered_env->getMinMz() - margin;
      if (min_mz < 0.0)  {
        min_mz = 0.0;
      }
      double max_mz = filtered_env->getMaxMz() + margin;
      of << feature->getId() << delimit
          << feature->getFracId() << delimit
          << single_feature->getEnvNum() << delimit
          << feature->getMonoMass() << delimit
          << mono_mz << delimit
          << charge << delimit
          << single_feature->getIntensity() << delimit
          << min_mz << delimit
          << max_mz << delimit
          << (single_feature->getTimeBegin()/60) << delimit
          << (single_feature->getTimeEnd()/60) << delimit
          << "#FF0000" << delimit
          << "0.1" << delimit
          << feature->getPromexScore()
          << std::endl;
    }
  }
  of.close();
}

void writeXmlFeatures(const std::string &output_file_name,
                      const FracFeaturePtrVec &features) {
  std::ofstream file;
  file.open(output_file_name.c_str());
  LOG_DEBUG("file_name " << output_file_name);
  file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  file << "<frac_feature_list>" << std::endl;

  for (size_t i = 0; i < features.size(); i++) {
    XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
    xercesc::DOMLSSerializer* serializer = impl->createSerializer();
    XmlDOMDocument doc(impl->createDoc("frac_feature_list"));
    XmlDOMElement* element = features[i]->toXmlElement(&doc);
    // LOG_DEBUG("Element generated");
    std::string str = xml_dom_util::writeToString(serializer, element);
    // LOG_DEBUG("String generated");
    xml_dom_util::writeToStreamByRemovingDoubleLF(file, str);
    element->release();
    serializer->release();
  }

  file << "</frac_feature_list>" << std::endl;
  file.close();
}

}

} /* namespace toppic */
