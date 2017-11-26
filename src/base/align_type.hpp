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


#ifndef PROT_BASE_ALIGN_TYPE_HPP_
#define PROT_BASE_ALIGN_TYPE_HPP_

#include <string>
#include <memory>

namespace prot {

class AlignType;
typedef std::shared_ptr<AlignType> AlignTypePtr;

class AlignType {
 public:
  static AlignTypePtr COMPLETE;
  static AlignTypePtr PREFIX;
  static AlignTypePtr SUFFIX;
  static AlignTypePtr INTERNAL;

  AlignType(const std::string &name, int id): name_(name), id_(id) {}

  std::string getName() {return name_;}

  int getId() {return id_;}

 private:
  std::string name_;
  int id_;
};

}  // namespace prot

#endif
