// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <boost/algorithm/string.hpp>

#include "base/logger.hpp"
#include "base/mass_constant.hpp"
#include "feature/envelope.hpp" 

namespace prot {

Envelope::Envelope(const Envelope &env) {
  refer_idx_ = env.refer_idx_;
  charge_ = env.charge_;
  mono_mz_ = env.mono_mz_;
  mzs_ = env.mzs_;
  intensities_ = env.intensities_;
}

Envelope::Envelope(int num, std::vector<std::string> &line_list) {
  charge_ = 1;
  std::vector<std::string> words;
  //LOG_DEBUG("line list size " << line_list.size() << " num " << num);
  boost::split(words, line_list[0], boost::is_any_of(" "));
  mono_mz_ = std::stod(words[7]);
  mzs_.resize(num);
  intensities_.resize(num);
  for (int i = 0; i < num; i++) {
    boost::split(words, line_list[i+1], boost::is_any_of(" "));
    mzs_[i] = std::stod(words[0]);
    intensities_[i] = std::stod(words[1])/100;
  }
  refer_idx_ = std::distance(intensities_.begin(), std::max_element(intensities_.begin(), intensities_.end()));
  //EnvelopeUtil::getMaxPos(intensities_);
}

Envelope::Envelope(int refer_idx, int charge, double mono_mz,
                   std::vector<double> &mzs, std::vector<double> &intensities):
    refer_idx_(refer_idx),
    charge_(charge),
    mono_mz_(mono_mz),
    mzs_(mzs),
    intensities_(intensities) {
    }

EnvelopePtr Envelope::convertToTheo(double mass_diff, int new_charge) {
  double new_mono_mz = (mono_mz_ + mass_diff) / new_charge;
  std::vector<double> new_mzs(mzs_.size());
  for (size_t i = 0; i < mzs_.size(); i++) {
    double new_value = (mzs_[i] + mass_diff)/ new_charge;
    new_mzs[i] = new_value;
  }
  EnvelopePtr env_ptr(new Envelope(refer_idx_, new_charge, new_mono_mz, 
                                   new_mzs, intensities_));
  return env_ptr;
}

// Convert a theoretical distribution to a theoretical envelope
EnvelopePtr Envelope::distrToTheoBase(double new_base_mz, int new_charge) {
  double mass_diff = new_base_mz * new_charge - mzs_[refer_idx_];
  return convertToTheo(mass_diff, new_charge);
}

// Convert a theoretical distribution to a theoretical envelope based on the
EnvelopePtr Envelope::distrToTheoMono(double new_mono_mz, int new_charge) {
  double mass_diff = new_mono_mz * new_charge - mono_mz_;
  return convertToTheo(mass_diff, new_charge);
}

void Envelope::changeIntensity(double ratio) {
  for (size_t i = 0; i < intensities_.size(); i++) {
    intensities_[i] *= ratio;
  }
}

void Envelope::changeToAbsInte(double absolute_intensity) {
  double ratio = absolute_intensity / getReferIntensity();
  changeIntensity(ratio);
}

void Envelope::changeMz(double shift) {
  for (size_t i = 0; i < mzs_.size(); i++) {
    mzs_[i] += shift;
  }
  mono_mz_ += shift;
}

EnvelopePtr Envelope::getSubEnv(int n_back, int n_forw) {
  int new_refer_idx = n_back;
  std::vector<double> new_mzs;
  std::vector<double> new_intes;
  for (int i = refer_idx_ - n_back; i <= refer_idx_ + n_forw; i++) {
    new_mzs.push_back(mzs_[i]);
    new_intes.push_back(intensities_[i]);
  }
  EnvelopePtr env_ptr(new Envelope(new_refer_idx, charge_, mono_mz_, 
                                   new_mzs, new_intes));
  return env_ptr;
}

EnvelopePtr Envelope::addZero(int num) {
  int n_peak = mzs_.size();
  std::vector<double> new_mzs(n_peak + 2 * num, 0); 
  std::vector<double> new_intes(n_peak + 2 * num, 0); 
  for (int i = 0; i < n_peak; i++) {
    new_mzs[i + num] = mzs_[i];
    new_intes[i + num] = intensities_[i];
  }
  for (int i = num - 1; i >= 0; i--) {
    new_mzs[i] = new_mzs[i+1] - MassConstant::getIsotopeMass() / charge_;
  }
  for (int i = n_peak + num; i < n_peak + num * 2; i++) {
    new_mzs[i] = new_mzs[i - 1] + MassConstant::getIsotopeMass() / charge_;
  }
  int new_refer_idx = refer_idx_ + num;
  EnvelopePtr env_ptr(new Envelope(new_refer_idx, charge_, mono_mz_, 
                                   new_mzs, new_intes));
  return env_ptr;
}

EnvelopePtr Envelope::getSubEnv(double percent_bound, double absolute_min_inte,
                                int max_back_peak_num, int max_forw_peak_num) {

  std::vector<int> bounds = calcBound(percent_bound, absolute_min_inte, 
                                      max_back_peak_num, max_forw_peak_num);
  return getSubEnv(bounds[0], bounds[1]);
}

// Compute the bound of highest peaks with intensity 85%.
std::vector<int> Envelope::calcBound(double percent_bound, double absolute_min_inte,
                                     int max_back_peak_num, int max_forw_peak_num) {
  double forw_inte;
  double back_inte;
  int forw_idx; // forward index 
  int back_idx; // backward index 
  // current sum of intensities
  double sum = intensities_[refer_idx_];
  forw_idx = refer_idx_ + 1;
  back_idx = refer_idx_ - 1;
  // if sum of intensity is above 85%, stop 
  double intensity_sum = compIntensitySum();
  int n_peak = mzs_.size();
  while (sum / intensity_sum < percent_bound) {
    // if all peaks are included, stop 
    if (back_idx == -1 && forw_idx == n_peak) {
      break;
    }
    // compute intensity at forward index 
    if (forw_idx == n_peak) {
      forw_inte = 0;
    } else {
      forw_inte = intensities_[forw_idx];
    }
    // compute intensity at backward index 
    if (back_idx == -1) {
      back_inte = 0;
    } else {
      back_inte = intensities_[back_idx];
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

void Envelope::shift(int shift) {
  refer_idx_ += shift;
  mono_mz_ += shift * MassConstant::getIsotopeMass() / charge_;
}

double Envelope::compIntensitySum() {
  double intensity_sum = 0;
  for (size_t i = 0; i < intensities_.size(); i++) {
    intensity_sum += intensities_[i];
  }
  return intensity_sum;
}

double Envelope::getAvgMz() {
  double sum = 0;
  for (size_t i = 0; i < mzs_.size(); i++) {
    if (mzs_[i] >= 0) {
      sum = sum + mzs_[i] * intensities_[i];
    }
  }
  return sum / compIntensitySum();
}

double Envelope::getAvgMass() {
  return getAvgMz() * charge_ - charge_ * MassConstant::getProtonMass();
}

}

