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


#ifndef PROT_SPEC_SPECTRUM_SET_HPP_
#define PROT_SPEC_SPECTRUM_SET_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_ms.hpp"
#include "spec/extend_ms.hpp"
#include "spec/prm_ms.hpp"
#include "spec/prm_ms_factory.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class SpectrumSet;
typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

class SpectrumSet {
 public:
  SpectrumSet(DeconvMsPtrVec deconv_ms_ptr_vec, SpParaPtr sp_para_ptr, 
              double prec_mono_mass);

  double getPrecMonoMass() {return prec_mono_mass_;}

  bool isValid() {return valid_;}

  int getSpectrumId() {return deconv_ms_ptr_vec_[0]->getMsHeaderPtr()->getId();}

  ExtendMsPtrVec getMsThreePtrVec() {return extend_ms_three_ptr_vec_;}

  DeconvMsPtrVec getDeconvMsPtrVec(){return deconv_ms_ptr_vec_;}

  PrmMsPtrVec getMsTwoPtrVec() {return prm_ms_two_ptr_vec_;}

  PrmMsPtrVec getSuffixMsTwoPtrVec() {return srm_ms_two_ptr_vec_;}

  PrmMsPtrVec getMsSixPtrVec(){return prm_ms_six_ptr_vec_;}

  PrmMsPtrVec getMsShiftSixPtrVec(double shift){
    return PrmMsFactory::geneShiftMsSixPtrVec(
        deconv_ms_ptr_vec_, sp_para_ptr_, prec_mono_mass_, -shift);
  }

 private:
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  SpParaPtr sp_para_ptr_;
  double prec_mono_mass_;
  bool valid_ = true;
  ExtendMsPtrVec extend_ms_three_ptr_vec_;
  PrmMsPtrVec prm_ms_two_ptr_vec_;
  PrmMsPtrVec srm_ms_two_ptr_vec_;
  PrmMsPtrVec prm_ms_six_ptr_vec_;

  bool checkValid(SpParaPtr sp_para_ptr);
};

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
