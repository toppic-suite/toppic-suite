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

#include "topfd/ecscore/spectrum/matrix_peak.hpp"

namespace toppic {

MatrixPeak::MatrixPeak(int peak_id, int spec_id, PeakPtr peak): 
  Peak(peak->getPosition(), peak->getIntensity()) {
    peak_id_ = peak_id;
    spec_id_ = spec_id;
    ori_inte_ = peak->getIntensity();
    start_idx_ = -1;
    end_idx_ = -1;
    neighbor_ = false;
  }

std::string MatrixPeak::getString() {
  return "Peak ID: " + std::to_string(peak_id_) + " " +
    "Spec ID: " + std::to_string(spec_id_) + " " +
    "Pos: " + std::to_string(getPosition()) + " " +
    "Inte: " + std::to_string(getIntensity()) + " " +
    "Neighbor: " + std::to_string(neighbor_) + "\n";
}

}
