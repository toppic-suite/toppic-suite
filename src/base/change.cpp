#include "change.hpp"

namespace prot {

Change::Change(int left_bp_pos, int right_bp_pos, int change_type,
         double mass_shift, PtmPtr ptm_ptr) {
  left_bp_pos_ = left_bp_pos;
  right_bp_pos_ = right_bp_pos;
  change_type_ = change_type;
  mass_shift_ = mass_shift;
  ptm_ptr_ = ptm_ptr;
}

}

