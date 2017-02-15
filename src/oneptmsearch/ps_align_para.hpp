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


#ifndef PROT_PS_ALIGN_PARA_HPP_
#define PROT_PS_ALIGN_PARA_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class PsAlignPara {
 public :
  PsAlignPara(int shift_num, double align_max_shift) { 
    n_unknown_shift_ = shift_num;
    align_max_shift_ = align_max_shift;
  }

  int getUnknownShiftNum() {return n_unknown_shift_;}

  /* parameters for ptm search */
  int n_unknown_shift_ =2;

  /* parameters for ps alignment */
  double align_max_shift_ = 1000000;
  // remove min shift so that large errors in precursor mass can be treated 
  // a as small shift
  //double align_min_shift_ = 0.5;
  
  // shift thresh for penalty
  double align_large_shift_thresh_ = 300;
  // set panelty to 3
  double align_large_shift_panelty_ = 3;

  double refine_prec_step_width_ = 0.005;

};

typedef std::shared_ptr<PsAlignPara> PsAlignParaPtr;

} /* namespace prot */

#endif /* PS_ALIGN_PARA_HPP_ */
