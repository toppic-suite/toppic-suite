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
#include <fstream>

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "ms/feature/frac_feature.hpp"

namespace toppic {

FracFeature::FracFeature(const std::string &file_name,
                         int frac_id, int feat_id,
                         double mono_mass, double inte,
                         int min_ms1_id, int max_ms1_id,
                         double time_begin, double time_end,
                         int scan_begin, int scan_end,
                         int min_charge, int max_charge, 
                         double apex_time, int apex_scan, 
                         double apex_inte, int rep_charge, 
                         double rep_avg_mz, int env_num,
                         double ec_score): 
    file_name_(file_name),
    frac_id_(frac_id),
    feat_id_(feat_id),
    mono_mass_(mono_mass),
    intensity_(inte),
    min_ms1_id_(min_ms1_id),
    max_ms1_id_(max_ms1_id),
    time_begin_(time_begin),
    time_end_(time_end),
    scan_begin_(scan_begin),
    scan_end_(scan_end),
    min_charge_(min_charge),
    max_charge_(max_charge),
    apex_time_(apex_time), 
    apex_scan_(apex_scan), 
    apex_inte_(apex_inte),
    rep_charge_(rep_charge),
    rep_avg_mz_(rep_avg_mz), 
    env_num_(env_num),
    ec_score_(ec_score) {
    }

FracFeature::FracFeature(std::string line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  file_name_ = strs[0];
  frac_id_ = std::stoi(strs[1]);
  feat_id_ = std::stoi(strs[2]);
  mono_mass_ = std::stod(strs[3]);
  intensity_ = std::stod(strs[4]);
  time_begin_ = std::stod(strs[5]) * 60;
  time_end_ = std::stod(strs[6]) * 60;
  scan_begin_ = std::stoi(strs[7]);
  scan_end_ = std::stoi(strs[8]);
  min_charge_ = std::stoi(strs[9]);
  max_charge_ = std::stoi(strs[10]);
  apex_time_ = std::stod(strs[11]) * 60;
  apex_scan_ = std::stoi(strs[12]);
  apex_inte_ = std::stod(strs[13]);
  rep_charge_ = std::stoi(strs[14]);
  rep_avg_mz_ = std::stod(strs[15]);
  env_num_ = std::stoi(strs[16]);
  ec_score_ = std::stod(strs[17]);
}

bool FracFeature::cmpFracIncInteDec(const FracFeaturePtr &a, 
                                    const FracFeaturePtr &b) { 
  if (a->getFracId() < b->getFracId()) {
    return true;
  }
  else if (a->getFracId() > b->getFracId()) {
    return false;
  }
  else if (a->getIntensity() > b->getIntensity()) {
    return true;
  }
  else {
    return false;
  }
}

FracFeature::FracFeature(XmlDOMElement* element) {
  file_name_ = xml_dom_util::getChildValue(element, "file_name", 0);
  frac_id_ = xml_dom_util::getIntChildValue(element, "frac_id", 0);
  feat_id_ = xml_dom_util::getIntChildValue(element, "feat_id", 0);
  mono_mass_ = xml_dom_util::getDoubleChildValue(element, "mono_mass", 0);
  intensity_ = xml_dom_util::getDoubleChildValue(element, "intensity", 0);
  time_begin_ = xml_dom_util::getDoubleChildValue(element, "time_begin", 0);
  time_end_ = xml_dom_util::getDoubleChildValue(element, "time_end", 0);
  scan_begin_ = xml_dom_util::getIntChildValue(element, "scan_begin", 0);
  scan_end_ = xml_dom_util::getIntChildValue(element, "scan_end", 0);
  min_charge_ = xml_dom_util::getIntChildValue(element, "min_charge", 0);
  max_charge_ = xml_dom_util::getIntChildValue(element, "max_charge", 0);
  apex_time_ = xml_dom_util::getDoubleChildValue(element, "apex_time", 0);
  apex_scan_ = xml_dom_util::getIntChildValue(element, "apex_scan", 0);
  apex_inte_ = xml_dom_util::getDoubleChildValue(element, "apex_inte", 0);
  rep_charge_ = xml_dom_util::getIntChildValue(element, "rep_charge", 0);
  rep_avg_mz_ = xml_dom_util::getDoubleChildValue(element, "rep_avg_mz", 0);
  env_num_ = xml_dom_util::getIntChildValue(element, "envelope_num", 0);
  ec_score_ = xml_dom_util::getDoubleChildValue(element, "ec_score", 0);

  // LOG_DEBUG("start parse changes");
  std::string single_feature_name = SingleChargeFeature::getXmlElementName();
  std::string feature_list_name = single_feature_name + "_list";

  XmlDOMElement* feature_list_element = xml_dom_util::getChildElement(element, feature_list_name.c_str(), 0);
  int feature_num = xml_dom_util::getChildCount(feature_list_element, single_feature_name.c_str());

  for (int i = 0; i < feature_num; i++) {
    XmlDOMElement* feature_element
        = xml_dom_util::getChildElement(feature_list_element, single_feature_name.c_str(), i);
    single_features_.push_back(std::make_shared<SingleChargeFeature>(feature_element));
  }
}


XmlDOMElement* FracFeature::toXmlElement(XmlDOMDocument* xml_doc) {
  std::string element_name = FracFeature::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "file_name", file_name_.c_str());
  std::string str = str_util::toString(frac_id_);
  xml_doc->addElement(element, "frac_id", str.c_str());
  str = str_util::toString(feat_id_);
  xml_doc->addElement(element, "feat_id", str.c_str());
  str = str_util::toString(mono_mass_);
  xml_doc->addElement(element, "mono_mass", str.c_str());
  str = str_util::toString(intensity_);
  xml_doc->addElement(element, "intensity", str.c_str());
  str = str_util::toString(time_begin_);
  xml_doc->addElement(element, "time_begin", str.c_str());
  str = str_util::toString(time_end_);
  xml_doc->addElement(element, "time_end", str.c_str());
  str = str_util::toString(scan_begin_);
  xml_doc->addElement(element, "scan_begin", str.c_str());
  str = str_util::toString(scan_end_);
  xml_doc->addElement(element, "scan_end", str.c_str());
  str = str_util::toString(min_charge_);
  xml_doc->addElement(element, "min_charge", str.c_str());
  str = str_util::toString(max_charge_);
  xml_doc->addElement(element, "max_charge", str.c_str());
  str = str_util::toString(apex_time_);
  xml_doc->addElement(element, "apex_time", str.c_str());
  str = str_util::toString(apex_scan_);
  xml_doc->addElement(element, "apex_scan", str.c_str());
  str = str_util::toString(apex_inte_);
  xml_doc->addElement(element, "apex_inte", str.c_str());
  str = str_util::toString(rep_charge_);
  xml_doc->addElement(element, "rep_charge", str.c_str());
  str = str_util::toString(rep_avg_mz_);
  xml_doc->addElement(element, "rep_avg_mz", str.c_str());
  str = str_util::toString(env_num_);
  xml_doc->addElement(element, "envelope_num", str.c_str());
  str = str_util::toString(ec_score_);
  xml_doc->addElement(element, "ec_score", str.c_str());

  element_name = SingleChargeFeature::getXmlElementName() + "_list";
  XmlDOMElement* cl = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < single_features_.size(); i++) {
    single_features_[i]->appendToXml(xml_doc, cl);
  }
  element->appendChild(cl);
  return element;
}


}
