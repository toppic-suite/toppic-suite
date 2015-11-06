#ifndef PROT_BASE_SEGMENT_HPP_
#define PROT_BASE_SEGMENT_HPP_

#include "base/trunc.hpp"

namespace prot {

class Segment {
 public:
  Segment(int left_bp_pos, int right_bp_pos, double n_shift, double c_shift);

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
