//
// Created by abbash on 8/26/22.
//

#include "simple_peak.hpp"

toppic::SimplePeak::SimplePeak(SimplePeak const &p) {
  pos_ = p.getPos();
  inte_ = p.getInte();
  start_idx_ = p.getStartIdx();
  end_idx_ = p.getEndIdx();
}

toppic::SimplePeak::SimplePeak(double pos, double inte) {
  pos_ = pos;
  inte_ = inte;
}