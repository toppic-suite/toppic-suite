// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <cmath>
#include <iostream>

#include "base/logger.hpp"
#include "feature/deconv_util.hpp"

namespace prot {

double DeconvUtil::intv_width_ = 10;

IntvDensPtrVec DeconvUtil::getDensity(std::vector<double> &inte) {
  double max_inte = -1;
  for (size_t i = 0; i < inte.size(); i++) {
    if (inte[i] > max_inte) {
      max_inte = inte[i];
    }
  }

  if (max_inte > 10000) {
    intv_width_ = max_inte / 1000;
  } else if (max_inte < 100) {
    intv_width_ = max_inte / 100;
  }
  int total_num = inte.size();
  int intv_num = (int)std::round(max_inte / intv_width_) + 1;
  IntvDensPtrVec dens;
  for (int i = 0; i < intv_num; i++) {
    double bgn = i * intv_width_;
    double end = (i + 1) * intv_width_;
    int num = 0;
    for (int j = 0; j < total_num; j++) {
      if (inte[j] > bgn && inte[j] <= end) {
        num++;
      }
    }
    IntvDensPtr cur_den(new IntvDens(bgn, end, num, num / (float) total_num));
    dens.push_back(cur_den);
  }
  return dens;
}

void DeconvUtil::outputDens(IntvDensPtrVec &dens) {
  for (size_t i = 0; i < dens.size(); i++) {
    std::cout << dens[i]->getBgn() << " " << dens[i]->getEnd() << " "
        << dens[i]->getNum() << " " << dens[i]->getPerc() << std::endl;
  }
}

int DeconvUtil::getMaxPos(IntvDensPtrVec &dens) {
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

double DeconvUtil::getBaseLine(std::vector<double> &inte) {
  //LOG_DEBUG("get density");
  IntvDensPtrVec dens = getDensity(inte);
  //LOG_DEBUG("get max pos ");
  int max_pos = getMaxPos(dens);
  LOG_DEBUG("max pos " << max_pos << " inte " << dens[max_pos]->getBgn());
  return dens[max_pos]->getBgn();
}

}
