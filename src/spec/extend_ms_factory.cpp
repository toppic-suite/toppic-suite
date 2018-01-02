//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <algorithm>
#include <vector>

#include "spec/extend_ms_factory.hpp"

namespace prot {

ExtendMsPtr ExtendMsFactory::geneMsThreePtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr,
                                            double new_prec_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mass);

  ExtendPeakPtrVec list;
  double ext_min_mass = sp_para_ptr->getExtendMinMass();
  std::vector<double> ext_offsets = sp_para_ptr->getExtendOffsets();
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    DeconvPeakPtr deconv_peak_ptr = deconv_ms_ptr->getPeakPtr(i);
    if (deconv_peak_ptr->getMonoMass() <= ext_min_mass) {
      double orig_mass = deconv_peak_ptr->getMonoMass();
      if (deconv_peak_ptr->getMonoMass() > new_prec_mass / 4) {
        orig_mass -= sp_para_ptr->mod_mass_[0];
      }

      if (deconv_peak_ptr->getMonoMass() > new_prec_mass * 3 / 4) {
        orig_mass -= sp_para_ptr->mod_mass_[1];
      }
      ExtendPeakPtr extend_peak_ptr
          = std::make_shared<ExtendPeak>(deconv_peak_ptr, orig_mass, 1.0);
      list.push_back(extend_peak_ptr);
    } else {
      for (size_t j = 0; j < ext_offsets.size(); j++) {
        double mass = deconv_peak_ptr->getMonoMass() + ext_offsets[j];
        if (deconv_peak_ptr->getMonoMass() + ext_offsets[j] > new_prec_mass / 4) {
          mass -= sp_para_ptr->mod_mass_[0];
        }

        if (deconv_peak_ptr->getMonoMass() + ext_offsets[j] > new_prec_mass * 3 / 4) {
          mass -= sp_para_ptr->mod_mass_[1];
        }
        ExtendPeakPtr extend_peak_ptr
            = std::make_shared<ExtendPeak>(deconv_peak_ptr, mass, 1.0);
        list.push_back(extend_peak_ptr);
      }
    }
  }

  // filter extend_peak
  ExtendPeakPtrVec list_filtered;
  double min_mass = sp_para_ptr->getMinMass();
  double prec_mono_mass = header_ptr->getPrecMonoMass();
  for (size_t i = 0; i < list.size(); i++) {
    double mass = list[i]->getPosition();
    if (mass >= min_mass && mass <= prec_mono_mass - min_mass) {
      list_filtered.push_back(list[i]);
    }
  }

  std::sort(list_filtered.begin(), list_filtered.end(), ExtendPeak::cmpPosIncrease);

  // set error tolerance
  PeakTolerancePtr peak_tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  for (size_t i = 0; i < list_filtered.size(); i++) {
    double mass = list_filtered[i]->getBasePeakPtr()->getMonoMass();
    double ori_tole = peak_tole_ptr->compStrictErrorTole(mass);
    list_filtered[i]->setOrigTolerance(ori_tole);
    double reve_tole = peak_tole_ptr->compRelaxErrorTole(mass, prec_mono_mass);
    list_filtered[i]->setReverseTolerance(reve_tole);
  }
  return std::make_shared<Ms<ExtendPeakPtr>>(header_ptr, list_filtered);
}

ExtendMsPtrVec ExtendMsFactory::geneMsThreePtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                                  SpParaPtr sp_para_ptr, double new_prec_mass) {
  ExtendMsPtrVec extend_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    extend_ms_ptr_vec.push_back(geneMsThreePtr(deconv_ms_ptr_vec[i], sp_para_ptr, new_prec_mass));
  }
  return extend_ms_ptr_vec;
}

} /* namespace prot */
