#ifndef PROT_CHANGE_HPP_
#define PROT_CHANGE_HPP_

#include "base/ptm.hpp"

namespace prot {

#define INPUT_CHANGE      0
#define FIXED_CHANGE      1
#define VARIABLE_CHANGE   2
#define UNEXPECTED_CHANGE 3

class Change {
 public:
  Change(int left_bp_pos, int right_bp_pos, int change_type,
         double mass_shift, PtmPtr ptm_ptr);
  
  int getLeftBpPos() {return left_bp_pos_;}

  int getRightBpPos() {return right_bp_pos_;}

  int getChangeType() {return change_type_;}

  double getMassShift() {return mass_shift_;}

  PtmPtr getPtmPtr() {return ptm_ptr_;}

 private:
  int left_bp_pos_;
  int right_bp_pos_;
  int change_type_;
  double mass_shift_;
  PtmPtr ptm_ptr_;
};

}
#endif

