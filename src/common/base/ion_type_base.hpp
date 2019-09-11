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


#ifndef TOPPIC_COMMON_BASE_ION_TYPE_BASE_HPP_
#define TOPPIC_COMMON_BASE_ION_TYPE_BASE_HPP_

#include "common/base/ion_type.hpp"

namespace toppic {

class IonTypeBase {
 public:
  static void initBase();

  static IonTypePtr getIonTypePtrByName(const std::string &name);

  static IonTypePtr getIonTypePtr_PREC() {return ion_type_ptr_PREC_;}

  static IonTypePtr getIonTypePtr_B() {return ion_type_ptr_B_;}

 private:
  static IonTypePtrVec ion_type_ptr_vec_;
  static IonTypePtr ion_type_ptr_B_;
  static IonTypePtr ion_type_ptr_PREC_;

  static std::string getName_B() {return "B";}
  static std::string getName_PREC() {return "PREC";}
};

}  // namespace toppic

#endif
