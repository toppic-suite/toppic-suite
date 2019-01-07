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

#include "deconv/feature/feature.hpp"

namespace toppic {

Feature::Feature(int id, double mono_mass, double inte,
                 int scan_begin, int scan_end): 
    id_(id),
    mono_mass_(mono_mass),
    intensity_(inte),
    scan_begin_(scan_begin),
    scan_end_(scan_end) {
    }


Feature::Feature(int id, double mono_mass, double inte,
                 double retent_begin, double retent_end,
                 int scan_begin, int scan_end,
                 int min_charge, int max_charge): 
    id_(id),
    mono_mass_(mono_mass),
    intensity_(inte),
    retent_begin_(retent_begin),
    retent_end_(retent_end),
    scan_begin_(scan_begin),
    scan_end_(scan_end),
    min_charge_(min_charge),
    max_charge_(max_charge) {
    }
}
