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


#ifndef PROT_FEATURE_REAL_ENVELOPE_HPP_
#define PROT_FEATURE_REAL_ENVELOPE_HPP_

#include <memory>
#include <vector>
#include <string>

#include "spec/peak.hpp"
#include "feature/envelope.hpp" 

namespace prot {

class RealEnv;

typedef std::shared_ptr<RealEnv> RealEnvPtr;

class RealEnv : public Envelope {
 public:
  RealEnv(std::vector<PeakPtr> &peak_list, EnvelopePtr theo_env, 
          double tolerance, double min_inte);

  int getNonExistPeakIdx() {return -1;}

  int getMissPeakNum() {return miss_peak_num_;}

  int getMatchPeakNum() {return getPeakNum() - miss_peak_num_;}

  int getMaxConsPeakNum() {return max_consecutive_peak_num_;}

  int getPeakIdx(int i) {return peak_idxes_[i];}

  std::vector<int> getPeakIdxList() {return peak_idxes_;}

  int getReferPeakIdx() {return peak_idxes_[refer_idx_];}

  static bool testPeakShare(RealEnvPtr a, RealEnvPtr b);

  bool isExist(int i);

 private:
  // peak index in the spectrum 
  // if peak_idx[i] == NO_EXIST_PEAK, it does not exist 
  std::vector<int> peak_idxes_;
  // number of missing peaks 
  int miss_peak_num_;
  // maximum number of consecutive peaks 
  int max_consecutive_peak_num_;

  void mapPeakList(std::vector<PeakPtr> &peak_list, EnvelopePtr theo_env, 
                   double tolerance, double min_min);

  void remvDuplMatch(EnvelopePtr theo_env);

  void cntMissPeakNum();

  void cntMaxConsPeakNum();
};

typedef std::vector<RealEnvPtr> RealEnvPtrVec;

}

#endif
