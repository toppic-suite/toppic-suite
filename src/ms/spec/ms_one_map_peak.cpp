//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#include "ms/spec/ms_one_map_peak.hpp"

namespace toppic {

MsOneMapPeak::MsOneMapPeak(double pos, double inte):
  Peak(pos, inte) {}


MsOneMapPeak::MsOneMapPeak(const MsOneMapPeakPtr p): 
  Peak(p->getPosition(), p->getIntensity()) { 
    start_idx_ = p->getStartIdx();
    end_idx_ = p->getEndIdx();
  }

}
