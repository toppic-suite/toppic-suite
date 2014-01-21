#ifndef PROT_SEGMENT_HPP_
#define PROT_SEGMENT_HPP_

#include "base/proteoform.hpp"
#include "base/trunc.hpp"

namespace prot {

class Segment {
 public:
  Segment(int left_bp_pos, int right_bp_pos, double n_shift, double c_shift) {
    left_bp_pos_ = left_bp_pos;
    right_bp_pos_ = right_bp_pos;
    pep_n_term_shift_ = n_shift;
    pep_c_term_shift_ = c_shift;
  }

  int getLeftBpPos () {return left_bp_pos_;}

  int getRightBpPos() {return right_bp_pos_;}

  double getPepNTermShift() {return pep_n_term_shift_;}

  double getPepCTermShift() {return pep_c_term_shift_;}

 private:
  // segment begin and end are based on break_points
  int left_bp_pos_;
  int right_bp_pos_;
  double pep_n_term_shift_; 
  double pep_c_term_shift_;
};

typedef std::shared_ptr<Segment> SegmentPtr;
typedef std::vector<SegmentPtr> SegmentPtrVec;

}
#endif
