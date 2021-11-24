//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "common/base/mass_constant.hpp"
#include "ms/spec/peak_util.hpp"

namespace toppic {

namespace peak_util {

double compPeakNeutralMass(double mono_mz, int charge) {
  return mono_mz * charge - charge * mass_constant::getProtonMass();
}

double compMz(double neutral_mass, int charge) {
    return neutral_mass / charge + mass_constant::getProtonMass();
}

}

}  // namespace toppic
