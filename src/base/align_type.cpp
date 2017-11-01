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

#include <string>

#include "base/align_type.hpp"

namespace prot {

AlignTypePtr AlignType::COMPLETE = AlignTypePtr(new AlignType("COMPLETE", 0));
AlignTypePtr AlignType::PREFIX = AlignTypePtr(new AlignType("PREFIX", 1));
AlignTypePtr AlignType::SUFFIX  = AlignTypePtr(new AlignType("SUFFIX", 2));
AlignTypePtr AlignType::INTERNAL = AlignTypePtr(new AlignType("INTERNAL", 3));

AlignType::AlignType(const std::string &name, int id):
    name_(name),
    id_(id) {
    }
}  // namespace prot
