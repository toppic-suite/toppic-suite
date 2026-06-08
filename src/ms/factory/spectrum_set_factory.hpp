//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_MS_FACTORY_SPECTRUM_SET_FACTORY_HPP_
#define TOPPIC_MS_FACTORY_SPECTRUM_SET_FACTORY_HPP_

#include <vector>

#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/spectrum_set.hpp"
#include "para/sp_para.hpp"

namespace toppic {

namespace spectrum_set_factory {

SpectrumSetPtr geneSpectrumSetPtr(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                  const SpParaPtr &sp_para_ptr,
                                  double prec_mono_mass);

SpectrumSetPtr readNextSpectrumSetPtr(const MsAlignReaderPtr &reader_ptr, 
                                      const SpParaPtr &sp_para_ptr);

SpectrumSetPtr readNextSpectrumSetPtr(const MsAlignReaderPtr &reader_ptr, 
                                      const SpParaPtr &sp_para_ptr, 
                                      int peak_num_limit);

SpectrumSetPtrVec geneSpectrumSetPtrVecWithPrecError(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                                     const SpParaPtr &sp_para_ptr,
                                                     const std::vector<double> &prec_error_vec);
}  // namespace spectrum_set_factory

}  // namespace toppic

#endif
