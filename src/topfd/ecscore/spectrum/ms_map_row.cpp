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

#include <iostream>

#include "topfd/ecscore/spectrum/ms_map_row.hpp"

namespace toppic {

MsMapRow::MsMapRow(MsMapRowHeaderPtr spec_ptr, int bin_num) {
    header_ptr_ = spec_ptr;
  for (int i = 0; i < bin_num; i++) {
    MsMapPeakPtrVec vec;
    peak_ptr_2d_.push_back(vec);
  }
}

void MsMapRow::print() {
  for (size_t i = 0; i < peak_ptr_2d_.size(); i++) {
    for (size_t j = 0; j < peak_ptr_2d_[i].size(); j++) {
      std::cout << "bin " << i << " j " << j << " m/z " << peak_ptr_2d_[i][j]->getPosition() << std::endl;
    }
  }
}
}
