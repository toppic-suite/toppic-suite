//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "spec/multi_spec_peak.hpp"

namespace prot {

MultiSpecPeak::MultiSpecPeak (double mono_mass, double intensity, 
                              int spec_id, int in_spec_peak_id, 
                              int base_type, double score):
    Peak (mono_mass, intensity),
    spec_id_(spec_id),
    in_spec_peak_id_(in_spec_peak_id),
    base_type_(base_type),
    score_(score) {
    }
}

