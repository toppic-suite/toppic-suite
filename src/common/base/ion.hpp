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

#ifndef TOPPIC_COMMON_BASE_ION_HPP_
#define TOPPIC_COMMON_BASE_ION_HPP_

#include "common/base/ion_type.hpp"
#include "common/base/neutral_loss.hpp"

namespace toppic {

class Ion {
 public:
  Ion(int charge, int pos, int display_pos,
      IonTypePtr ion_type_ptr,
      NeutralLossPtr neutral_loss_ptr);

  int getCharge() const {return charge_;}

  int getPos() const {return pos_;}

  int getDisplayPos() const {return display_pos_;}

  IonTypePtr getIonTypePtr() const {return ion_type_ptr_;}

  std::string getDisplayName() const;

 private:
  int charge_;
  int pos_;
  int display_pos_;
  IonTypePtr ion_type_ptr_;
  NeutralLossPtr neutral_loss_ptr_;
};

using IonPtr = std::shared_ptr<Ion>;
using IonPtrVec = std::vector<IonPtr>;

}  // namespace toppic

#endif 
