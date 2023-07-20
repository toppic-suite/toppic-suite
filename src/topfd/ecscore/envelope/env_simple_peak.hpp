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

#ifndef TOPPIC_TOPFD_ECSCORE_ENVELOPE_ENV_SIMPLE_PEAK_HPP
#define TOPPIC_TOPFD_ECSCORE_ENVELOPE_ENV_SIMPLE_PEAK_HPP

#include "ms/spec/peak.hpp"

namespace toppic {

class EnvSimplePeak;

typedef std::shared_ptr<EnvSimplePeak> EnvSimplePeakPtr;

class EnvSimplePeak : public Peak {
 public:
  EnvSimplePeak();

  EnvSimplePeak(double pos, double inte);

  EnvSimplePeak(const EnvSimplePeakPtr p);

  int getStartIdx() const { return start_idx_; }

  void setStartIdx(int startIdx) { start_idx_ = startIdx; }

  int getEndIdx() const { return end_idx_; }

  void setEndIdx(int endIdx) { end_idx_ = endIdx; }

  static bool cmpInteDec(EnvSimplePeakPtr a, EnvSimplePeakPtr b) {
    return a->getIntensity() > b->getIntensity(); }

  static bool cmpPosInc(EnvSimplePeakPtr a, EnvSimplePeakPtr b) {
    return a->getPosition() < b->getPosition(); }

 private:
  int start_idx_ = -1;
  int end_idx_ = -1;
};

typedef std::vector<EnvSimplePeakPtr> EnvSimplePeakPtrVec;

}

#endif 
