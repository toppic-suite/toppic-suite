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


#ifndef TOPPIC_COMMON_BASE_SUPPORT_PEAK_TYPE_BASE_HPP_
#define TOPPIC_COMMON_BASE_SUPPORT_PEAK_TYPE_BASE_HPP_

#include "common/base/support_peak_type.hpp"

namespace toppic {

class SPTypeBase {
 public:
  static void initBase();

  static const SPTypePtrVec& getBaseSPTypePtrVec() {
    return sp_type_ptr_vec_;}

  static SPTypePtr getSPTypePtrByName(const std::string &name);
  static SPTypePtr getSPTypePtrById(int id);

  static SPTypePtr getSPTypePtr_N_TERM() {
    return sp_type_ptr_N_TERM_;
  }

  static std::string getName_N_TERM() {return "N_TERM";}

 private:
  static SPTypePtrVec sp_type_ptr_vec_;
  static SPTypePtr sp_type_ptr_N_TERM_;
};

} /* namespace toppic */

#endif /* SUPPORT_PEAK_TYPE_BASE_HPP_ */
