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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/env/env_base.hpp"
#include "topfd/ecscore/env/seed_env.hpp"

namespace toppic {

SeedEnv::SeedEnv(DeconvPeakPtr peak_ptr) {
  //init for parent class
  mono_mz_ = peak_ptr->getMonoMz();
  charge_ = peak_ptr->getCharge();
  double mass = peak_ptr->getMonoMass();
  EnvPtr theo_env_ptr = EnvBase::getEnvByMonoMass(mass, charge_);
  refer_idx_ = theo_env_ptr->getReferIdx();
  for (int i = 0; i < theo_env_ptr->getPeakNum(); i++) {
    EnvPeakPtr p_ptr = std::make_shared<EnvPeak>(theo_env_ptr->getMz(i),
                                                 theo_env_ptr->getInte(i));
    peak_ptr_list_.push_back(p_ptr);
  }
  // init for seed envelope
  spec_id_ = peak_ptr->getSpId();
  seed_inte_ = peak_ptr->getIntensity();
}

SeedEnv::SeedEnv(SeedEnvPtr env_ptr, int new_charge) {
  //init for parent class
  charge_ = new_charge;
  double mass = env_ptr->getMonoNeutralMass(); 
  mono_mz_ = peak_util::compMz(mass, charge_);
  refer_idx_ = env_ptr->getReferIdx();
  int old_charge = env_ptr->getCharge();
  for (auto &i: env_ptr->peak_ptr_list_) {
    double old_mz = i->getPosition();
    double neutral_mass = peak_util::compPeakNeutralMass(old_mz, old_charge);
    double new_mz = peak_util::compMz(neutral_mass, new_charge); 
    EnvPeakPtr p_ptr = std::make_shared<EnvPeak>(new_mz,
                                                 i->getIntensity());
    peak_ptr_list_.push_back(p_ptr);
  }
  // init for seed envelope
  spec_id_ = env_ptr->spec_id_;
  seed_inte_ = env_ptr->seed_inte_;
}

SeedEnv::SeedEnv(SeedEnvPtr env_ptr, EnvPeakPtrVec &peak_ptr_list) {
  //init for parent class
  mono_mz_ = env_ptr->getMonoMz(); 
  charge_ = env_ptr->getCharge();
  int left_side_non_exist_num = 0;
  for (size_t i = 0; i < peak_ptr_list.size(); i++) {
    if (peak_ptr_list[i]->isExist()) {
      break;
    }
    else {
      left_side_non_exist_num++;
    }
  }
  refer_idx_ = env_ptr->getReferIdx() - left_side_non_exist_num;
  for (size_t i = 0; i < peak_ptr_list.size(); i++) {
    if (peak_ptr_list[i]->isExist()) {
      peak_ptr_list_.push_back(peak_ptr_list[i]);
    }
  }
  // init for seed envelope
  spec_id_ = env_ptr->spec_id_;
  seed_inte_ = env_ptr->seed_inte_;
}

SeedEnv::SeedEnv(SeedEnvPtr env_ptr) {
  //init for parent class
  mono_mz_ = env_ptr->getMonoMz(); 
  refer_idx_ = env_ptr->getReferIdx();
  charge_ = env_ptr->getCharge();
  for (auto &i: env_ptr->peak_ptr_list_) {
    EnvPeakPtr p_ptr = std::make_shared<EnvPeak>(i->getPosition(), 
                                                 i->getIntensity());
    peak_ptr_list_.push_back(p_ptr);
  }
  // init for seed envelope
  spec_id_ = env_ptr->spec_id_;
  seed_inte_ = env_ptr->seed_inte_;
}


EnvPeakPtrVec SeedEnv::getScaledPeakPtrList(double ratio, double min_inte) {
  EnvPeakPtrVec new_peak_list; 
  int non_exist_idx = EnvPeak::getNonExistPeakIdx();
  double non_exist_mz = -1; 
  double non_exist_inte = 0;
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    EnvPeakPtr peak_ptr = peak_ptr_list_[i];
    double scaled_inte = peak_ptr->getIntensity() * ratio;
    if (scaled_inte < min_inte) {
      EnvPeakPtr new_peak_ptr = std::make_shared<EnvPeak>(non_exist_mz,
                                                          non_exist_inte,
                                                          non_exist_idx);
      new_peak_list.push_back(new_peak_ptr);
    }
    else {
      EnvPeakPtr new_peak_ptr = std::make_shared<EnvPeak>(peak_ptr->getPosition(), 
                                                          scaled_inte, i);
      new_peak_list.push_back(new_peak_ptr);
    }
  }
  return new_peak_list;
}

std::string SeedEnv::getString() {
  std::string header = "Spec ID: " + std::to_string(spec_id_) + " " +
                       "Pos: " + std::to_string(mono_mz_) + " " +
                       "Inte: " + std::to_string(seed_inte_) + " " +
                       "Charge: " + std::to_string(charge_) + "\n";
  std::string peaks = "(";
  for (auto peak: peak_ptr_list_)
    peaks = peaks + "(" + std::to_string(peak->getPosition()) + ", " 
      + std::to_string(peak->getIntensity()) + "), ";
  peaks = peaks + ")\n";
  return header + peaks;
}

void SeedEnv::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = "seed_envelope";
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(spec_id_);
  xml_doc->addElement(element, "spec_id", str.c_str());
  str = str_util::toString(seed_inte_);
  xml_doc->addElement(element, "seed_inte", str.c_str());
  str = str_util::toString(getMonoNeutralMass());
  xml_doc->addElement(element, "mono_mass", str.c_str());
  str = str_util::toString(getCharge());
  xml_doc->addElement(element, "charge", str.c_str());
  str = str_util::toString(getReferMz());
  xml_doc->addElement(element, "refer_mz", str.c_str());
  element_name = "peak_list";
  XmlDOMElement* peak_list = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    peak_ptr_list_[i]->appendToXml(xml_doc, peak_list);
  }
  element->appendChild(peak_list);
  parent->appendChild(element);
}

}
