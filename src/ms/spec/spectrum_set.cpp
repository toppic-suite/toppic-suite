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
#include "ms/spec/extend_ms_factory.hpp"
#include "ms/spec/prm_ms_factory.hpp"
#include "ms/spec/spectrum_set.hpp"

namespace toppic {

SpectrumSet::SpectrumSet(DeconvMsPtrVec deconv_ms_ptr_vec, 
                         double prec_mono_mass,
                         bool valid, 
                         ExtendMsPtrVec extend_ms_three_ptr_vec,
                         PrmMsPtrVec prm_ms_two_ptr_vec,
                         PrmMsPtrVec srm_ms_two_ptr_vec,
                         PrmMsPtrVec prm_ms_six_ptr_vec): 
    deconv_ms_ptr_vec_(deconv_ms_ptr_vec),
    prec_mono_mass_(prec_mono_mass), 
    valid_(valid),
    extend_ms_three_ptr_vec_(extend_ms_three_ptr_vec),
    prm_ms_two_ptr_vec_(prm_ms_two_ptr_vec),
    srm_ms_two_ptr_vec_(prm_ms_two_ptr_vec),
    prm_ms_six_ptr_vec_(prm_ms_six_ptr_vec) {}

}  // namespace toppic
