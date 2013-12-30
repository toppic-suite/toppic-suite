#include "base/change.hpp"

namespace prot {

Change::Change(int left_bp_pos, int right_bp_pos, int change_type,
         double mass_shift, PtmPtr ptm_ptr) {
  left_bp_pos_ = left_bp_pos;
  right_bp_pos_ = right_bp_pos;
  change_type_ = change_type;
  mass_shift_ = mass_shift;
  ptm_ptr_ = ptm_ptr;
}

Change::Change(Change ori, int start) {
  left_bp_pos_ = ori.left_bp_pos_ - start;
  right_bp_pos_ = ori.right_bp_pos_ - start;
  change_type_ = ori.change_type_;
  mass_shift_ = ori.mass_shift_;
  ptm_ptr_ = ori.ptm_ptr_;
}

}

