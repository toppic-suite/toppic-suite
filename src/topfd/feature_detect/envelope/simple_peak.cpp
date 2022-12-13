//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#include "simple_peak.hpp"

namespace toppic {
  SimplePeak::SimplePeak() {
    pos_ = -1;
    inte_ = -1;
  }

  SimplePeak::SimplePeak(SimplePeak const &p) {
    pos_ = p.getPos();
    inte_ = p.getInte();
    start_idx_ = p.getStartIdx();
    end_idx_ = p.getEndIdx();
  }

  SimplePeak::SimplePeak(double pos, double inte) {
    pos_ = pos;
    inte_ = inte;
  }

  bool SimplePeak::isEmpty() {
    if (pos_ == -1 and inte_ == -1) return true;
    return false;
  }
}