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


#include <cmath>
#include <iostream>
#include <vector>

#include "base/logger.hpp"
#include "feature/deconv_util.hpp"

namespace prot {

namespace deconv_util {

IntvDensPtrVec getDensity(const std::vector<double> &inte) {
  double max_inte = -1;
  for (size_t i = 0; i < inte.size(); i++) {
    if (inte[i] > max_inte) {
      max_inte = inte[i];
    }
  }

  double intv_width = 10;

  if (max_inte > 10000) {
    intv_width = max_inte / 1000;
  } else if (max_inte < 100) {
    intv_width = max_inte / 100;
  }
  size_t total_num = inte.size();
  int intv_num = static_cast<int>(std::round(max_inte / intv_width)) + 1;
  IntvDensPtrVec dens(intv_num);
  for (int i = 0; i < intv_num; i++) {
    double bgn = i * intv_width;
    double end = (i + 1) * intv_width;
    int num = 0;
    for (size_t j = 0; j < total_num; j++) {
      if (inte[j] > bgn && inte[j] <= end) {
        num++;
      }
    }
    IntvDensPtr cur_den
        = std::make_shared<IntvDens>(bgn, end, num, num / static_cast<float>(total_num));
    dens[i] = cur_den;
  }
  return dens;
}

void outputDens(const IntvDensPtrVec &dens) {
  for (size_t i = 0; i < dens.size(); i++) {
    std::cout << dens[i]->getBgn() << " " << dens[i]->getEnd() << " "
        << dens[i]->getNum() << " " << dens[i]->getPerc() << std::endl;
  }
}

int getMaxPos(const IntvDensPtrVec &dens) {
  int max_pos = -1;
  int max_num = -1;
  for (size_t i = 0; i < dens.size(); i++) {
    if (dens[i]->getNum() > max_num) {
      max_num = dens[i]->getNum();
      max_pos = i;
    }
  }
  return max_pos;
}

double getBaseLine(const std::vector<double> &inte) {
  // LOG_DEBUG("get density");
  IntvDensPtrVec dens = getDensity(inte);
  // LOG_DEBUG("get max pos ");
  int max_pos = getMaxPos(dens);
  LOG_DEBUG("max pos " << max_pos << " inte " << dens[max_pos]->getBgn());
  return dens[max_pos]->getBgn();
}

} // namespace deconv_util

}  // namespace prot
