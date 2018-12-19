//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_BASE_SEGMENT_HPP_
#define TOPPIC_BASE_SEGMENT_HPP_

#include <vector>

#include "base/trunc.hpp"

namespace toppic {

class Segment {
 public:
  Segment(int left_bp_pos, int right_bp_pos, double n_shift, double c_shift):
      left_bp_pos_(left_bp_pos),
      right_bp_pos_(right_bp_pos),
      pep_n_term_shift_(n_shift),
      pep_c_term_shift_(c_shift) {}

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

}  // namespace toppic
#endif
