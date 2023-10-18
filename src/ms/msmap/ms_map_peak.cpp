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

#include "ms/msmap/ms_map_peak.hpp"

namespace toppic {

MsMapPeak::MsMapPeak(PeakPtr peak):
  Peak(peak->getPosition(), peak->getIntensity()) {
    ori_inte_ = peak->getIntensity();
    neighbor_ = false;
  }

std::string MsMapPeak::getString() {
  return "Pos: " + std::to_string(getPosition()) + " " +
    "Inte: " + std::to_string(getIntensity()) + " " +
    "Neighbor: " + std::to_string(neighbor_) + "\n";
}

}
