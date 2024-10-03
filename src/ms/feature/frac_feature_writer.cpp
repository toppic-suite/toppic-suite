// Copyright (c) 2014 - 2023, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <algorithm>
#include <iomanip>
#include <cstddef>
#include <set>
#include <sstream>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_impl.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "ms/env/env.hpp"
#include "ms/env/env_base.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/feature/frac_feature_writer.hpp"

namespace toppic {

namespace frac_feature_writer {

void writeHeader(std::ofstream &of) {
  of.precision(16);
  of << "File_name" << "\t"
     << "Fraction_ID" << "\t"
     << "Feature_ID" << "\t"
     << "Mass" << "\t"
     << "Intensity" << "\t"
     << "Min_time" << "\t"
     << "Max_time" << "\t"
     << "Min_scan" << "\t"
     << "Max_scan" << "\t"
     << "Min_charge" << "\t"
     << "Max_charge" << "\t"
     << "Apex_time" << "\t"
     << "Apex_scan" << "\t"
     << "Apex_intensity" << "\t"
     << "Rep_charge" << "\t"
     << "Rep_average_mz" << "\t"
     << "Envelope_num" << "\t"
     << "EC_score" << "\t" << std::endl;
}

void writeOneFeature(std::ofstream &of, FracFeaturePtr feature) {
  of << feature->getFileName() << "\t" << feature->getFracId() << "\t"
     << feature->getFeatId() << "\t" << feature->getMonoMass() << "\t"
     << feature->getIntensity() << "\t" << feature->getTimeBegin() / 60 << "\t"
     << feature->getTimeEnd() / 60 << "\t" << feature->getScanBegin() << "\t"
     << feature->getScanEnd() << "\t" << feature->getMinCharge() << "\t"
     << feature->getMaxCharge() << "\t" << feature->getApexTime() / 60 << "\t"
     << feature->getApexScan() << "\t" << feature->getApexInte() << "\t"
     << feature->getRepCharge() << "\t" << feature->getRepAvgMz() << "\t"
     << feature->getEnvNum() << "\t" << feature->getEcScore() << "\t"
     << std::endl;
}

void writeFeatures(const std::string &output_file_name,
                   const FracFeaturePtrVec &features) {
  std::ofstream of(output_file_name);
  writeHeader(of);

  for (std::size_t i = 0; i < features.size(); i++) {
    FracFeaturePtr feature = features[i];
    writeOneFeature(of, feature);
  }
  of.close();
}

void writeBatMassFeatures(const std::string &output_file_name,
                          const FracFeaturePtrVec &features, int num_spectra) {
  std::ofstream of(output_file_name);
  std::string delimit = ",";
  of << "ID" << delimit << "Fraction_ID" << delimit << "Envelope_num" << delimit
     << "Mass" << delimit << "MonoMz" << delimit << "Charge" << delimit
     << "Intensity" << delimit << "mzLo" << delimit << "mzHi" << delimit
     << "rtLo" << delimit << "rtHi" << delimit << "specLow" << delimit
     << "specHi" << delimit << "rtApex" << delimit << "Score" << delimit
     << "XIC" << delimit << "Envelope" << std::endl;
  for (std::size_t i = 0; i < features.size(); i++) {
    FracFeaturePtr feature = features[i];
    double mono_mass = feature->getMonoMass();
    SingleChargeFeaturePtrVec single_features = feature->getSingleFeatures();

    for (std::size_t j = 0; j < single_features.size(); j++) {
      SingleChargeFeaturePtr single_feature = single_features[j];
      int charge = single_feature->getCharge();
      double mono_mz = peak_util::compMz(mono_mass, charge);
      EnvPtr ref_env = EnvBase::getEnvByMonoMass(mono_mass);
      EnvPtr theo_env = ref_env->distrToTheoMono(mono_mz, charge);
      double min_inte = 0.03;
      EnvPtr filtered_env = theo_env->getSubEnv(min_inte);
      /////
      std::vector<double> xic = single_feature->getXicInte();
      std::stringstream ss;
      int k = 0;
      for (int i = 0; i < num_spectra; i++) {
        if (i != 0) ss << ";";
        if (i >= single_feature->getSpecIDBegin() and
            i <= single_feature->getSpecIDEnd()) {
          ss << xic[k];
          k++;
        } else
          ss << 0;
      }

      std::vector<double> envelopeMass = single_feature->getEnvelopeMass();
      std::vector<double> aggregateEnvelopeInte =
          single_feature->getAggregateEnvelopeInte();
      std::stringstream env;
      for (std::size_t i = 0; i < envelopeMass.size(); i++) {
        if (i != 0) env << ";";
        env << envelopeMass[i];
        env << '&';
        env << aggregateEnvelopeInte[i];
      }
      /////
      // margin for envelopes
      double margin = 0.1;
      double min_mz = filtered_env->getMinMz() - margin;
      if (min_mz < 0.0) {
        min_mz = 0.0;
      }
      double max_mz = filtered_env->getMaxMz() + margin;
      of << feature->getFeatId() << delimit << feature->getFracId() << delimit
         << single_feature->getEnvNum() << delimit;
      of << std::fixed << std::setprecision(5) << feature->getMonoMass()
         << delimit << mono_mz << delimit << charge << delimit
         << single_feature->getIntensity() << delimit << min_mz << delimit
         << max_mz << delimit;
      of << std::setprecision(2) << (single_feature->getTimeBegin() / 60)
         << delimit << (single_feature->getTimeEnd() / 60) << delimit
         << single_feature->getSpecIDBegin() << delimit
         << single_feature->getSpecIDEnd() << delimit
         << feature->getApexTime() / 60 << delimit << feature->getEcScore()
         << delimit << ss.str() << delimit << env.str() << std::endl;
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

  for (std::size_t i = 0; i < features.size(); i++) {
    XmlDOMImpl *impl = XmlDOMImplFactory::getXmlDOMImplInstance();
    xercesc::DOMLSSerializer *serializer = impl->createSerializer();
    XmlDOMDocument doc(impl->createDoc("frac_feature_list"));
    XmlDOMElement *element = features[i]->toXmlElement(&doc);
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

}  // namespace frac_feature_writer

} /* namespace toppic */
