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
#include "ms/spec/deconv_ms_util.hpp"
#include "ms/spec/extend_ms_factory.hpp"
#include "ms/spec/prm_ms_factory.hpp"
#include "ms/spec/spectrum_set_factory.hpp"

namespace toppic {

namespace spectrum_set_factory {

bool checkValid(DeconvMsPtrVec &deconv_ms_ptr_vec, SpParaPtr sp_para_ptr,
                double prec_mono_mass) {
  if (prec_mono_mass < sp_para_ptr->getMinMass()) {
    return false;
  }
  int peak_num = 0;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    peak_num += deconv_ms_ptr_vec[i]->size();
  }

  LOG_DEBUG("peak_num " << peak_num);
  if (peak_num < sp_para_ptr->getMinPeakNum()) {
    return false;
  }
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    if (deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr() == nullptr) {
      return false;
    }
  }
  return true;
}

SpectrumSetPtr geneSpectrumSetPtr(DeconvMsPtrVec deconv_ms_ptr_vec,
                                  SpParaPtr sp_para_ptr,
                                  double prec_mono_mass) { 
  bool valid = checkValid(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);

  ExtendMsPtrVec extend_ms_three_ptr_vec;

  PrmMsPtrVec prm_ms_two_ptr_vec;

  PrmMsPtrVec srm_ms_two_ptr_vec;

  PrmMsPtrVec prm_ms_six_ptr_vec;

  if (valid) {
    extend_ms_three_ptr_vec
        = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
    prm_ms_two_ptr_vec
        = prm_ms_factory::geneMsTwoPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
    srm_ms_two_ptr_vec
        = prm_ms_factory::geneSuffixMsTwoPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
    prm_ms_six_ptr_vec
        = prm_ms_factory::geneMsSixPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
  }
  SpectrumSetPtr spec_set_ptr = std::make_shared<SpectrumSet>(deconv_ms_ptr_vec, 
                                                              prec_mono_mass, 
                                                              valid, 
                                                              extend_ms_three_ptr_vec,
                                                              prm_ms_two_ptr_vec,
                                                              srm_ms_two_ptr_vec,
                                                              prm_ms_six_ptr_vec);
  return spec_set_ptr;
}

SpectrumSetPtr readNextSpectrumSetPtr(SimpleMsAlignReaderPtr reader_ptr, 
                                      SpParaPtr sp_para_ptr) {
  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  if (deconv_ms_ptr_vec.size() == 0) {
    return nullptr;
  }
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();

  SpectrumSetPtr spec_set_ptr = geneSpectrumSetPtr(deconv_ms_ptr_vec,
                                                   sp_para_ptr, prec_mono_mass);
  return spec_set_ptr;
}

SpectrumSetPtr readNextSpectrumSetPtr(SimpleMsAlignReaderPtr reader_ptr, 
                                      SpParaPtr sp_para_ptr, 
                                      int peak_num_limit) {
  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  if (deconv_ms_ptr_vec.size() == 0) {
    return nullptr;
  }
  // keep only top peaks
  deconv_ms_util::keepTopPeaks(deconv_ms_ptr_vec, peak_num_limit);
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  SpectrumSetPtr spec_set_ptr = geneSpectrumSetPtr(deconv_ms_ptr_vec,
                                                   sp_para_ptr, prec_mono_mass);
  return spec_set_ptr;
}


SpectrumSetPtrVec geneSpectrumSetPtrVecWithPrecError(DeconvMsPtrVec deconv_ms_ptr_vec,  
                                                     SpParaPtr sp_para_ptr) {
  SpectrumSetPtrVec spec_set_vec;
  if (deconv_ms_ptr_vec.size() == 0) {
    return spec_set_vec;
  }
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  std::vector<double> prec_error_vec = sp_para_ptr->getPrecErrorVec(); 
  for (size_t i = 0; i< prec_error_vec.size(); i++) {
    SpectrumSetPtr spec_set_ptr = geneSpectrumSetPtr(deconv_ms_ptr_vec,
                                                     sp_para_ptr, 
                                                     prec_mono_mass + prec_error_vec[i]);
    spec_set_vec.push_back(spec_set_ptr);
  }
  return spec_set_vec;
}

}
}  // namespace toppic
