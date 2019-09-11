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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "feature/frac_xml_feature_reader.hpp"

namespace toppic {

FracXmlFeatureReader::FracXmlFeatureReader(const std::string &file_name):
    file_name_(file_name) {
      input_.open(file_name.c_str(), std::ios::in);
      if (!input_.is_open()) {
        LOG_ERROR("Feature file  " << file_name << " does not exist.");
        exit(EXIT_FAILURE); 
      }
      // read header line
      std::string line;
      std::getline(input_, line);
    }

FracXmlFeatureReader::~FracXmlFeatureReader() {
  if (input_.is_open()) {
    input_.close();
  }
}

void FracXmlFeatureReader::close() {
  input_.close();
}

std::vector<std::string> FracXmlFeatureReader::readOneFeatureLines() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    str_util::trim(line);
    // LOG_DEBUG("line " << line);
    if (line ==  "<frac_feature>") {
      line_list.push_back(line);
    } else if (line == "</frac_feature>") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } else if (line == "") {
      continue;
    } else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

FracFeaturePtr FracXmlFeatureReader::readOneFeature() {
  std::vector<std::string> ft_str_vec = readOneFeatureLines();
  if (ft_str_vec.size() == 0) {
    return FracFeaturePtr(nullptr);
  }
  std::string ft_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
  for (size_t i = 0; i < ft_str_vec.size(); i++) {
    ft_str += ft_str_vec[i];
  }
  // LOG_DEBUG("prsm str " << prsm_str);
  xercesc::MemBufInputSource ft_buf(
      (const XMLByte*)ft_str.c_str(), ft_str.size(), "feature_str (in memory)");

  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  FracFeaturePtr ptr;
  if (parser) {
    XmlDOMDocument doc(parser, ft_buf);
    XmlDOMElement* root = doc.getDocumentElement();
    ptr = std::make_shared<FracFeature>(root);
  }
  return ptr;
}

FracFeaturePtrVec FracXmlFeatureReader::readAllFeatures() {
  FracFeaturePtrVec all_features;
  FracFeaturePtr feature;
  while ((feature = readOneFeature()) != nullptr) {
    all_features.push_back(feature);
  }
  return all_features;
}


}  // namespace prot
