//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_MS_SPEC_BASE_PEAK_TYPE_HPP_
#define TOPPIC_MS_SPEC_BASE_PEAK_TYPE_HPP_

#include <memory>
#include <string>

namespace toppic {

class BasePeakType;
typedef std::shared_ptr<BasePeakType> BasePeakTypePtr;

class BasePeakType {
 public:
  static const BasePeakTypePtr ORIGINAL;
  static const BasePeakTypePtr REVERSED;

  BasePeakType(std::string name) {name_=name;}

  std::string getName() {return name_;}

 private:
  std::string name_;
};

} // namespace toppic

#endif

