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


#include <string>
#include <vector>

#include "util/str_util.hpp"
#include "base/ion.hpp"

namespace toppic {

Ion::Ion(int charge, int pos, int display_pos,
         IonTypePtr ion_type_ptr,
         NeutralLossPtr neutral_loss_ptr):
    charge_(charge),
    pos_(pos),
    display_pos_(display_pos),
    ion_type_ptr_(ion_type_ptr),
    neutral_loss_ptr_(neutral_loss_ptr) {}

std::string Ion::getDisplayName() {
  return ion_type_ptr_->getName() + str_util::toString(display_pos_);
}

} /* namespace toppic */

