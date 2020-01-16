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


#ifndef TOPPIC_SEARCH_ONE_PTM_SEARCH_PS_ALIGN_PARA_HPP_
#define TOPPIC_SEARCH_ONE_PTM_SEARCH_PS_ALIGN_PARA_HPP_

#include "prsm/prsm_para.hpp"

namespace toppic {

class PsAlignPara {
 public :
  PsAlignPara(int shift_num, double align_max_shift, double align_min_shift): 
      n_unknown_shift_(shift_num),
      align_max_shift_(align_max_shift),
      align_min_shift_(align_min_shift) { }

  int getUnknownShiftNum() {return n_unknown_shift_;}

  /* parameters for ptm search */
  int n_unknown_shift_ =2;

  /* parameters for ps alignment */
  double align_max_shift_ = 1000000;

  double align_min_shift_ = -1000000;

  // shift thresh for penalty
  double align_large_shift_thresh_ = 300;
  // set panelty to 3
  double align_large_shift_panelty_ = 3;

  double refine_prec_step_width_ = 0.005;

};

typedef std::shared_ptr<PsAlignPara> PsAlignParaPtr;

} /* namespace toppic */

#endif /* PS_ALIGN_PARA_HPP_ */
