// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "base/logger.hpp"
#include "feature/raw_ms_util.hpp"
#include "feature/deconv_data_base.hpp"

namespace prot {

DeconvDataPtr DeconvDataBase::getDataPtr(PeakPtrVec &peak_list, FeatureMngPtr mng_ptr) {
  if (peak_list.size() == 0) return nullptr;
  double max_mz = RawMsUtil::findMaxPos(peak_list);
  if (max_mz > mng_ptr->max_mass_) {
    LOG_WARN("Max mz is too large: " << max_mz);
    return nullptr;
  }
  return std::make_shared<DeconvData>(peak_list, mng_ptr->max_mass_,
                                      mng_ptr->max_charge_, mng_ptr->window_size_);
}


// generate deconvolution data using given max mass, max charge
DeconvDataPtr DeconvDataBase::getDataPtr(PeakPtrVec &peak_list, double max_mass,
                                         int max_charge, FeatureMngPtr mng_ptr) {
  if (max_charge < 1) {
    LOG_WARN("Max charge < 1");
    max_charge = mng_ptr->max_charge_;
  }
  if (max_mass <= 0) {
    LOG_WARN("Max mass <= 0");
    max_mass = mng_ptr->max_mass_;
  }
  if (max_mass > mng_ptr->max_mass_) {
    LOG_WARN("Max mass is greater than default max mass " << max_mass);
    max_mass = mng_ptr->max_mass_;
  }
  double max_mz = RawMsUtil::findMaxPos(peak_list);
  if (max_mz > mng_ptr->max_mass_) {
    LOG_WARN("Max mz is too large: " << max_mz);
    return nullptr;
  }
  for (size_t i = 0; i < peak_list.size(); i++) {
    if (peak_list[i]->getPosition() < 0 || peak_list[i]->getIntensity() < 0) {
      LOG_WARN("mz intensity are negative values");
      return nullptr;
    }
  }
  return std::make_shared<DeconvData>(peak_list, max_mass,
                                      max_charge, mng_ptr->window_size_);
}

}  // namespace prot
