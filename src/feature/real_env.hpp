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
