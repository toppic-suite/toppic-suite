#ifndef PROT_CHANGE_HPP_
#define PROT_CHANGE_HPP_

#include "base/ptm.hpp"

namespace prot {

#define INPUT_CHANGE      0
#define FIXED_CHANGE      1
#define PROTEIN_VARIABLE_CHANGE   2
#define VARIABLE_CHANGE   3
#define UNEXPECTED_CHANGE 4

class Change {
 public:
  Change(int left_bp_pos, int right_bp_pos, int change_type,
         double mass_shift, PtmPtr ptm_ptr);

  Change(Change ori_change, int shift);
  
  int getLeftBpPos() {return left_bp_pos_;}

  int getRightBpPos() {return right_bp_pos_;}

  int getChangeType() {return change_type_;}

  double getMassShift() {return mass_shift_;}

  PtmPtr getPtmPtr() {return ptm_ptr_;}

 private:
  // left and right positions are based on break point positions 
  int left_bp_pos_;
  int right_bp_pos_;
  int change_type_;
  double mass_shift_;
  PtmPtr ptm_ptr_;
};

typedef std::shared_ptr<Change> ChangePtr;
typedef std::vector<ChangePtr> ChangePtrVec;

inline bool compareChangeUp(ChangePtr c1, ChangePtr c2) {
  return c1->getLeftBpPos() < c2->getLeftBpPos();
}

}
#endif

