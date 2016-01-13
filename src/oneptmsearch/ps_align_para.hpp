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
