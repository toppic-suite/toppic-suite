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


#ifndef TOPPIC_SEQ_BREAK_POINT_HPP_
#define TOPPIC_SEQ_BREAK_POINT_HPP_

#include <memory>
#include <vector>

#include "base/ion_type.hpp"

namespace toppic {

class BreakPoint {
 public:
  BreakPoint(double prm, double srm);

  double getPrm() {return prm_;}

  double getSrm() {return srm_;}

  double getNTermMass(IonTypePtr ion_type_ptr);

  double getCTermMass(IonTypePtr ion_type_ptr);

 private:
  double prm_;
  double srm_;
};

typedef std::shared_ptr<BreakPoint> BreakPointPtr;
typedef std::vector<BreakPointPtr> BreakPointPtrVec;

}  // namespace toppic

#endif
