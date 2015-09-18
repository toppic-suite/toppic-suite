#ifndef PROT_LOCAL_PROCESSOR_HPP_
#define PROT_LOCAL_PROCESSOR_HPP_

#include "htslib/faidx.h"
#include "base/ptm.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "spec/theo_peak.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "local_mng.hpp"

namespace prot {

class LocalProcessor {

  public:
    LocalProcessor(LocalMngPtr& mng_ptr);
    void process();
    void init();

  private:
    void processOneSpectrum(PrsmPtr &prsm, PrsmWriter & writer);

    void processOnePtm(PrsmPtr& prsm);
    void processOneKnown(PrsmPtr& prsm);
    void processOneUnknown(PrsmPtr& prsm);

    void processTwoPtm(PrsmPtr& prsm);
    void processTwoKnown(PrsmPtr& prsm);
    void processTwoUnknown(PrsmPtr& prsm, bool ptm1_known, bool ptm2_known);

    LocalMngPtr mng_ptr_;
    faidx_t *fai_;
    std::vector<double> para_;
    double p1_, p2_, ppm_;
    double weight_, theta_, thread_, beta_;
    double max_ptm_mass_, min_mass_;
    bool cysteine_protected_;
    PrsmPtrVec prsm_ptrs_;
};

typedef std::shared_ptr<LocalProcessor> LocalProcessorPtr;


}

#endif
