//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include "common/base/mod_type.hpp"

namespace toppic {

ModTypePtr ModType::SIDE_CHAIN = std::make_shared<ModType>("SIDE_CHAIN", 0);

ModTypePtr ModType::N_TERM = std::make_shared<ModType>("N_TERM", 1);

ModTypePtr ModType::C_TERM = std::make_shared<ModType>("C_TERM", 2);

ModType::ModType(const std::string &name, int id): 
    name_(name), 
    id_(id) {}

ModTypePtr ModType::getModTypePtrByName(std::string name) {
  if (name == "SIDE_CHAIN") {
    return ModType::SIDE_CHAIN;
  }
  if (name == "N_TERM") {
    return ModType::N_TERM;
  }
  if (name == "C_TERM") {
    return ModType::C_TERM;
  }
  return nullptr;
}

}  // namespace toppic
