//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <cmath>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/base/mass_constant.hpp"
#include "ms/env/env.hpp"

namespace toppic {

Env::Env(Env &env):
    refer_idx_(env.refer_idx_),
    charge_(env.charge_),
    mono_mz_(env.mono_mz_) {
      for (int i = 0; i < env.getPeakNum(); i++) {
        EnvPeakPtr peak_ptr = std::make_shared<EnvPeak>(env.getPeakPtr(i));
        peak_ptr_list_.push_back(peak_ptr);
      }
    }

Env::Env(int num, std::vector<std::string> &line_list) {
  charge_ = 1;
  std::vector<std::string> words = str_util::split(line_list[0], " ");
  peak_ptr_list_.resize(num);
  for (int i = 0; i < num; i++) {
    words = str_util::split(line_list[i+1], " ");
    double mz = std::stod(words[0]);
    double inte = std::stod(words[1]) / 100;
    EnvPeakPtr peak_ptr = std::make_shared<EnvPeak>(mz, inte);
      peak_ptr_list_[i] = peak_ptr;
  }
  refer_idx_ = getHighestPeakIdx();
  mono_mz_ = peak_ptr_list_[0]->getPosition();
}

Env::Env(int refer_idx, int charge, double mono_mz,
         EnvPeakPtrVec &peaks):
    refer_idx_(refer_idx),
    charge_(charge),
    mono_mz_(mono_mz) {
      for (size_t i = 0; i < peaks.size(); i++) {
        EnvPeakPtr peak_ptr = std::make_shared<EnvPeak>(peaks[i]);
        peak_ptr_list_.push_back(peak_ptr);
      }
    }

EnvPtr Env::convertToTheo(double mass_diff, int new_charge) {
  int ori_charge = 1;
  double new_mono_mz = (mono_mz_ * ori_charge + mass_diff) / new_charge;
  EnvPeakPtrVec new_peaks(peak_ptr_list_.size());
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    double new_mz = (peak_ptr_list_[i]->getPosition() + mass_diff) / new_charge;
    new_peaks[i] = std::make_shared<EnvPeak>(new_mz, peak_ptr_list_[i]->getIntensity());
  }
  return std::make_shared<Env>(refer_idx_, new_charge, new_mono_mz,
                               new_peaks);
}

// Convert a theoretical distribution to a theoretical envelope
EnvPtr Env::distrToTheoRef(double new_ref_mz, int new_charge) {
  int ori_charge = 1;
  double mass_diff = new_ref_mz * new_charge - peak_ptr_list_[refer_idx_]->getPosition() * ori_charge;
  return convertToTheo(mass_diff, new_charge);
}

// Convert a theoretical distribution to a theoretical envelope based on the
EnvPtr Env::distrToTheoMono(double new_mono_mz, int new_charge) {
  int ori_charge = 1;
  double mass_diff = new_mono_mz * new_charge - mono_mz_ * ori_charge;
  return convertToTheo(mass_diff, new_charge);
}

void Env::changeIntensity(double ratio) {
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    peak_ptr_list_[i]->changeIntensity(ratio);
  }
}

void Env::changeToAbsInte(double absolute_intensity) {
  double ratio = absolute_intensity / getReferInte();
  changeIntensity(ratio);
}

void Env::changeMz(double shift) {
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    peak_ptr_list_[i]->shiftPosition(shift);
  }
  mono_mz_ += shift;
}

void Env::changeMzByIsotope(double shift_num) {
  double shift_mass = shift_num * mass_constant::getIsotopeMass();
  double shift_mz = shift_mass / charge_;
  changeMz(shift_mz);
}

EnvPtr Env::getSubEnv(int n_back, int n_forw) {
  int new_refer_idx = n_back;
  EnvPeakPtrVec new_peaks;
  for (int i = refer_idx_ - n_back; i <= refer_idx_ + n_forw; i++) {
    new_peaks.push_back(peak_ptr_list_[i]);
  }
  EnvPtr env_ptr = std::make_shared<Env>(new_refer_idx, charge_, mono_mz_,
                                         new_peaks);
  return env_ptr;
}

EnvPtr Env::addZero(int num) {
  int n_peak = peak_ptr_list_.size();
  EnvPeakPtrVec new_peaks; 
  for (int i = 0; i < n_peak + 2 * num; i++) {
    EnvPeakPtr peak_ptr = std::make_shared<EnvPeak>(0 , 0);
    new_peaks.push_back(peak_ptr);
  }
  for (int i = 0; i < n_peak; i++) {
    new_peaks[i + num]->setPosition(peak_ptr_list_[i]->getPosition());
    new_peaks[i + num]->setIntensity(peak_ptr_list_[i]->getIntensity());
  }
  for (int i = num - 1; i >= 0; i--) {
    double pos = new_peaks[i+1]->getPosition() 
        - mass_constant::getIsotopeMass() / charge_;
    new_peaks[i]->setPosition(pos);
  }
  for (int i = n_peak + num; i < n_peak + num * 2; i++) {
    double pos = new_peaks[i-1]->getPosition() 
        + mass_constant::getIsotopeMass() / charge_;
    new_peaks[i]->setPosition(pos);
  }
  int new_refer_idx = refer_idx_ + num;
  EnvPtr env_ptr = std::make_shared<Env>(new_refer_idx, charge_, mono_mz_,
                                         new_peaks);
  return env_ptr;
}

EnvPtr Env::getSubEnv(double percent_bound, double absolute_min_inte,
                      int max_back_peak_num, int max_forw_peak_num) {
  std::vector<int> bounds = calcBound(percent_bound, absolute_min_inte,
                                      max_back_peak_num, max_forw_peak_num);
  return getSubEnv(bounds[0], bounds[1]);
}

EnvPtr Env::getSubEnv(double min_inte) {
  size_t left = refer_idx_;
  size_t right = refer_idx_;
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    if (peak_ptr_list_[i]->getIntensity() >= min_inte) {
      if (i < left) {
        left = i;
      }
      if (i > right) {
        right = i;
      }
    }
  }
  return getSubEnv(refer_idx_ - left, right - refer_idx_); 
}


// Compute the bound of highest peaks with intensity 85%.
std::vector<int> Env::calcBound(double percent_bound, double absolute_min_inte,
                                int max_back_peak_num, int max_forw_peak_num) {
  double forw_inte;
  double back_inte;
  int forw_idx;  // forward index
  int back_idx;  // backward index
  // current sum of intensities
  double sum = peak_ptr_list_[refer_idx_]->getIntensity();
  forw_idx = refer_idx_ + 1;
  back_idx = refer_idx_ - 1;
  // if sum of intensity is above 85%, stop
  double intensity_sum = compInteSum();
  int n_peak = peak_ptr_list_.size();
  while (sum / intensity_sum < percent_bound) {
    // if all peaks are included, stop
    if (back_idx == -1 && forw_idx == n_peak) {
      break;
    }
    // compute intensity at forward index
    if (forw_idx == n_peak) {
      forw_inte = 0;
    } else {
      forw_inte = peak_ptr_list_[forw_idx]->getIntensity();
    }
    // compute intensity at backward index
    if (back_idx == -1) {
      back_inte = 0;
    } else {
      back_inte = peak_ptr_list_[back_idx]->getIntensity();
    }
    // if both forw and back intensities are less than a threshold, stop
    if (back_inte < absolute_min_inte && forw_inte < absolute_min_inte) {
      break;
    }
    if (back_inte > forw_inte) {
      sum += back_inte;
      back_idx--;
    } else {
      sum += forw_inte;
      forw_idx++;
    }
  }
  int max_back_num = refer_idx_ - back_idx - 1;
  int max_forw_num = forw_idx - refer_idx_ - 1;

  if (max_back_num > max_back_peak_num) {
    max_back_num = max_back_peak_num;
  }
  if (max_forw_num > max_forw_peak_num) {
    max_forw_num = max_forw_peak_num;
  }
  // return results
  std::vector<int> result;
  result.push_back(max_back_num);
  result.push_back(max_forw_num);
  return result;
}

void Env::changeReferIdx(int shift) {
  refer_idx_ += shift;
  mono_mz_ += shift * mass_constant::getIsotopeMass() / charge_;
}

double Env::compInteSum() {
  double sum = 0;
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    sum = sum + peak_ptr_list_[i]->getIntensity();
  }
  return sum;
}

double Env::getAvgMz() {
  double sum = 0;
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    if (peak_ptr_list_[i]->getPosition() >= 0) {
      sum = sum + peak_ptr_list_[i]->getPosition() * peak_ptr_list_[i]->getIntensity();
    }
  }
  return sum / compInteSum();
}

int Env::getHighestPeakIdx() {
  return std::distance(peak_ptr_list_.begin(),
                       std::max_element(peak_ptr_list_.begin(), peak_ptr_list_.end(), EnvPeak::cmpInteInc));
}

std::vector<double> Env::getInteList() {
  std::vector<double> intensities;
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    intensities.push_back(peak_ptr_list_[i]->getIntensity());
  }
  return intensities;
}

std::vector<double> Env::getMzList() {
  std::vector<double> pos_list;
  for (auto p: peak_ptr_list_)
    pos_list.push_back(p->getPosition());
  return pos_list;
}


void Env::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = Env::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(refer_idx_);
  xml_doc->addElement(element, "refer_idx", str.c_str());
  str = str_util::toString(charge_);
  xml_doc->addElement(element, "charge", str.c_str());
  str = str_util::toString(mono_mz_);
  xml_doc->addElement(element, "mono_mz", str.c_str());
  for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
    peak_ptr_list_[i]->appendXml(xml_doc, element);
  }
  parent->appendChild(element);
}

}  // namespace toppic
