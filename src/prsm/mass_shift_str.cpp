//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "prsm/mass_shift_str.hpp"

namespace toppic {

MassShiftStr::MassShiftStr(double mass_shift, int left_pos, 
                           int right_pos):
    mass_shift_(mass_shift),
    left_pos_(left_pos),
    right_pos_(right_pos) {}


bool MassShiftStr::cmpPosInc(const std::shared_ptr<MassShiftStr> &a,
                             const std::shared_ptr<MassShiftStr> &b) {
  if (a->left_pos_ < b->left_pos_) {
    return true;
  } else if (a->left_pos_ > b->left_pos_) {
    return false;
  } else {
    return a->right_pos_ < b->right_pos_;
  }
}


}  // namespace toppic

