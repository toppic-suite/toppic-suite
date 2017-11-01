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


#ifndef PROT_BASE_ACTIVATION_BASE_HPP_
#define PROT_BASE_ACTIVATION_BASE_HPP_

#include <string>

#include "base/activation.hpp"

namespace prot {

class ActivationBase {
 private:
  static ActivationPtrVec activation_ptr_vec_;

 public:
  static void initBase(const std::string &file_name);

  static const ActivationPtrVec& getActivationPtrVec() {return activation_ptr_vec_;}

  static ActivationPtr getActivationPtrByName(const std::string &name);

  static ActivationPtr getActivationPtrFromXml(xercesc::DOMElement * element);
};

}  // namespace prot

#endif
