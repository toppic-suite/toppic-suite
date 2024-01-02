//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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


#include "common/util/logger.hpp"
#include "topfd/dia/isolation_window.hpp"

namespace toppic {

IsolationWindow::IsolationWindow(double mz_bgn, double mz_end):
  mz_bgn_(mz_bgn),
  mz_end_(mz_end) {}

bool IsolationWindow::findSpecId(int id) {
  if (spec_id_set_.find(id) == spec_id_set_.end()) {
    return false;
  }
  else {
    return true;
  }
}

bool IsolationWindow::isMatch(double mz_bgn, double mz_end) {
  if (mz_bgn == mz_bgn_ && mz_end == mz_end_) {
    return true;
  }
  else {
    return false;
  }
}

}  // namespace toppic
