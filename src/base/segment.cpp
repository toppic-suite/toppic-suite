#include "base/segment.hpp"

namespace prot {

Segment::Segment(int left_bp_pos, int right_bp_pos, 
                 double n_shift, double c_shift):
    left_bp_pos_(left_bp_pos),
    right_bp_pos_(right_bp_pos),
    pep_n_term_shift_(n_shift),
    pep_c_term_shift_(c_shift) {
    }

}
