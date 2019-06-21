//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_DECONV_ENV_REAL_ENVELOPE_HPP_
#define TOPPIC_DECONV_ENV_REAL_ENVELOPE_HPP_

#include "spec/peak.hpp"
#include "deconv/env/envelope.hpp" 

namespace toppic {

class RealEnv;

typedef std::shared_ptr<RealEnv> RealEnvPtr;

class RealEnv : public Envelope {
 public:
  RealEnv(const PeakPtrVec &peak_list, EnvelopePtr theo_env, 
          double tolerance, double min_inte);

  int getSpId() {return sp_id_;}

  int getMissPeakNum() {return miss_peak_num_;}

  int getMatchPeakNum() {return getPeakNum() - miss_peak_num_;}

  int getMaxConsPeakNum() {return max_consecutive_peak_num_;}

  int getPeakIdx(int i) {return peaks_[i]->getIdx();}

  int getReferPeakIdx() {return peaks_[refer_idx_]->getIdx();}

  void setSpId(int sp_id) {sp_id_ = sp_id;}

  static bool testPeakShare(RealEnvPtr a, RealEnvPtr b);

  bool isExist(int i);

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "real_env";}

 private:
  int sp_id_;
  // peak index in the spectrum 
  // if peak_idx[i] == NO_EXIST_PEAK, it does not exist 
  // std::vector<int> peak_idxes_;
  // number of missing peaks 
  int miss_peak_num_;
  // maximum number of consecutive peaks 
  int max_consecutive_peak_num_;

  void mapPeakList(const PeakPtrVec &peak_list, EnvelopePtr theo_env, 
                   double tolerance, double min_min);

  void remvDuplMatch(EnvelopePtr theo_env);

  void cntMissPeakNum();

  void cntMaxConsPeakNum();
};

typedef std::vector<RealEnvPtr> RealEnvPtrVec;
typedef std::vector<RealEnvPtrVec> RealEnvPtrVec2D;

}

#endif
