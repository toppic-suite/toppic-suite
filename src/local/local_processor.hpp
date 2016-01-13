#ifndef PROT_LOCAL_PROCESSOR_HPP_
#define PROT_LOCAL_PROCESSOR_HPP_

#include "htslib/faidx.h"
#include "base/ptm.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "spec/theo_peak.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "local_mng.hpp"
#include "local_anno.hpp"

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
  double theta_, thread_, beta_;
  double min_mass_;
};

typedef std::shared_ptr<LocalProcessor> LocalProcessorPtr;

}
#endif
