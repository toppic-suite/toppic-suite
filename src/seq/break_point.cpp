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


#include "common/base/ion_type.hpp"
#include "seq/break_point.hpp"

namespace toppic {

BreakPoint::BreakPoint(double prm, double srm): 
    prm_(prm), srm_(srm) {}

double BreakPoint::getNTermMass(IonTypePtr ion_type_ptr) {
  return prm_ + ion_type_ptr->getShift();
}

double BreakPoint::getCTermMass(IonTypePtr ion_type_ptr) {
  return srm_ + ion_type_ptr->getShift();
}

}  // namespace toppic
