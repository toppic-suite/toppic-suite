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

#include "common/util/logger.hpp"
#include "para/peak_tolerance_base.hpp"

namespace toppic {

PeakTolerancePtrVec PeakToleranceBase::peak_tole_ptr_vec_;

PeakTolerancePtr PeakToleranceBase::main_peak_tole_ptr_;

std::unordered_map<std::string, PeakTolerancePtr> PeakToleranceBase::peak_tole_name_map_;

void PeakToleranceBase::initBase(std::string name, double ppm) {
  double ppo = ppm *0.000001;
  bool use_min_tolerance = true;
  double min_tolerance = 0.01;
  main_peak_tole_ptr_ 
      = std::make_shared<PeakTolerance>(name, ppo, use_min_tolerance, min_tolerance);

  peak_tole_ptr_vec_.push_back(main_peak_tole_ptr_);

  peak_tole_name_map_[main_peak_tole_ptr_->getName()] = main_peak_tole_ptr_;
}


PeakTolerancePtr PeakToleranceBase::getPeakTolerancePtrByName(const std::string &name) {
  PeakTolerancePtr peak_tole_ptr = peak_tole_name_map_[name];
  if (peak_tole_ptr == nullptr) {
    LOG_ERROR("Peak tolerance " << name << " cannot be found!");
  }
  return peak_tole_ptr;
}

}  // namespace toppic
