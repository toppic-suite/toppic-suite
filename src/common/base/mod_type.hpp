//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_COMMON_BASE_MOD_TYPE_HPP_
#define TOPPIC_COMMON_BASE_MOD_TYPE_HPP_

#include <string>
#include <memory>
#include <vector>

namespace toppic {

class ModType;
using ModTypePtr = std::shared_ptr<ModType>;

class ModType {
 public:
  static ModTypePtr SIDE_CHAIN;
  static ModTypePtr N_TERM;
  static ModTypePtr C_TERM;

  ModType(const std::string &name, int id);

  const std::string& getName() const {return name_;}

  int getId() const {return id_;}

  static ModTypePtr getModTypePtrByName(const std::string &name);

 private:
  std::string name_;
  int id_;
};

using ModTypePtrVec = std::vector<ModTypePtr>;

}  // namespace toppic

#endif
