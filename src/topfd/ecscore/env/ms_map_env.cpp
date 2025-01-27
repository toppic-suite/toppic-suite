//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include "topfd/ecscore/env/ms_map_env.hpp"

namespace toppic {

MsMapEnv::MsMapEnv(int spec_id, MsMapPeakPtrVec peak_list) {
  spec_id_ = spec_id;
  peak_list_ = peak_list;
}

int MsMapEnv::getTopThreeMatchNum(int ref_idx) {
  int total_peaks = peak_list_.size();
  int start_idx = std::max(ref_idx - 1, 0);
  int end_idx = std::min(ref_idx + 1, total_peaks - 1);
  int num = 0;
  for (int i = start_idx; i < end_idx + 1; i++) {
    MsMapPeakPtr p = peak_list_[i];
    if (p != nullptr) {
      num = num + 1;
    }
  }
  return num;
}

std::vector<double> MsMapEnv::getInteList() {
  std::vector<double> inte_list;
  for (auto p: peak_list_) {
    if (p != nullptr)
      inte_list.push_back(p->getIntensity());
    else
      inte_list.push_back(0);
  }
  return inte_list;
}

std::vector<double> MsMapEnv::getMzList() {
  std::vector<double> pos_list;
  for (auto p: peak_list_) {
    if (p != nullptr)
      pos_list.push_back(p->getPosition());
    else
      pos_list.push_back(0);
  }
  return pos_list;
}

double MsMapEnv::getInteSum() {
  double sum = 0;
  for (auto p: peak_list_) {
    if (p != nullptr)
      sum = sum + p->getIntensity();
  }
  return sum;
}


double MsMapEnv::compTopThreeInteSum(int ref_idx) {
  double sum = 0;
  int peak_num = peak_list_.size();
  if (ref_idx >= 0 && ref_idx < peak_num 
      && peak_list_[ref_idx] != nullptr) {
    sum += peak_list_[ref_idx]->getIntensity();
  }
  int left_idx = ref_idx - 1;
  if (left_idx >= 0 && left_idx < peak_num 
      && peak_list_[left_idx] != nullptr) {
    sum += peak_list_[left_idx]->getIntensity();
  }
  int right_idx = ref_idx + 1;
  if (right_idx >= 0 && right_idx < peak_num 
      && peak_list_[right_idx] != nullptr) {
    sum += peak_list_[right_idx]->getIntensity();
  }
  return sum;
}

void MsMapEnv::removeLowIntePeaks(SeedEnvPtr seed_ptr, double ratio,
                                  double min_inte) {
  for (int i = 0; i < seed_ptr->getPeakNum(); i++) {
    double inte = seed_ptr->getPeakPtr(i)->getIntensity() * ratio;
    if (inte < min_inte) {
      peak_list_[i] = nullptr;
    }
  }
}

void MsMapEnv::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = "ms_map_envelope";
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(spec_id_);
  xml_doc->addElement(element, "spec_id", str.c_str());
  str = str_util::toString(getInteSum());
  xml_doc->addElement(element, "inte_sum", str.c_str());
  element_name = "peak_list";
  XmlDOMElement* peak_list_elem = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < peak_list_.size(); i++) {
    if (peak_list_[i] != nullptr) {
      peak_list_[i]->appendToXml(xml_doc, peak_list_elem);
    }
  }
  element->appendChild(peak_list_elem);
  parent->appendChild(element);
}

}

