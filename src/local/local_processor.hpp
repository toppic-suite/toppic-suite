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


#ifndef PROT_LOCAL_PROCESSOR_HPP_
#define PROT_LOCAL_PROCESSOR_HPP_

#include "htslib/faidx.h"
#include "base/ptm.hpp"
#include "base/local_anno.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "spec/theo_peak.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "local_mng.hpp"

namespace prot {

class LocalProcessor {

 public:
  LocalProcessor(LocalMngPtr mng_ptr);
  void process();

 private:
  void processOneSpectrum(PrsmPtr prsm);

  void processOnePtm(PrsmPtr prsm);
  ProteoformPtr processOneKnown(const PrsmPtr & prsm);
  ProteoformPtr processOneUnknown(const PrsmPtr & prsm);

  void processTwoPtm(PrsmPtr prsm);
  ProteoformPtr processTwoKnown(const PrsmPtr & prsm);
  ProteoformPtr processTwoUnknown(const PrsmPtr & prsm);

  LocalMngPtr mng_ptr_;
  double p1_, p2_, ppm_;
  double threshold_; // threshold for MIScore;
  double theta_; // the weight for known/unknown ptm
  double beta_; // the weight for one/two ptm
  double min_mass_;
};

typedef std::shared_ptr<LocalProcessor> LocalProcessorPtr;

}
#endif
