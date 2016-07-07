#include <cmath>
#include <iostream>

#include "base/logger.hpp"
#include "feature/deconv_util.hpp"

namespace prot {

double intv_width_ = 10;

IntvDensPtrVec getDensity(std::vector<double> &inte) {
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
  LOG_DEBUG("get density");
  IntvDensPtrVec dens = getDensity(inte);
  LOG_DEBUG("get max pos ");
  int max_pos = getMaxPos(dens);
  return dens[max_pos]->getEnd();
}

}
