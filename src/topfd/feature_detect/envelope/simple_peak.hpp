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

#ifndef TOPPIC_SIMPLE_PEAK_HPP
#define TOPPIC_SIMPLE_PEAK_HPP

#include "topfd/feature_detect/spectrum/peak_matrix.hpp"

namespace toppic {
  class SimplePeak {
  public:
    SimplePeak();

    SimplePeak(double pos, double inte);

    SimplePeak(SimplePeak const &p);

    bool isEmpty();

    double getPos() const { return pos_; }

    void setPos(double pos) { pos_ = pos; }

    double getInte() const { return inte_; }

    void setInte(double inte) { inte_ = inte; }

    int getStartIdx() const { return start_idx_; }

    void setStartIdx(int startIdx) { start_idx_ = startIdx; }

    int getEndIdx() const { return end_idx_; }

    void setEndIdx(int endIdx) { end_idx_ = endIdx; }

    static bool cmpInteDec(const SimplePeak &a, const SimplePeak &b) { return a.getInte() > b.getInte(); }

    static bool cmpPos(const SimplePeak &a, const SimplePeak &b) { return a.getPos() < b.getPos(); }

  private:
    double pos_;
    double inte_;
    int start_idx_ = -1;
    int end_idx_ = -1;

  };
}

#endif //TOPPIC_SIMPLE_PEAK_HPP
