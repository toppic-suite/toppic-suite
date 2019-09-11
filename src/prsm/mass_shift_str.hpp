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

#ifndef TOPPIC_PRSM_MASS_SHIFT_STR_HPP_
#define TOPPIC_PRSM_MASS_SHIFT_STR_HPP_

#include <memory>

namespace toppic {

class MassShiftStr {
 public:
  MassShiftStr(double mass_shift, int left_pos, int right_pos);

  double mass_shift_;
  int left_pos_; 
  int right_pos_;
  static bool cmpPosInc(const std::shared_ptr<MassShiftStr> &a,
                        const std::shared_ptr<MassShiftStr> &b);
};

typedef std::shared_ptr<MassShiftStr> MassShiftStrPtr;

}  // namespace toppic

#endif
