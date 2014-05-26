#ifndef PROT_POSSION_PROCESSOR_HPP_
#define PROT_POSSION_PROCESSOR_HPP_

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "poisson/poisson_mng.hpp"
#include "poisson/poisson_comp_pvalue.hpp"

namespace prot {

class PoissonProcessor {
 public:
  PoissonProcessor(PoissonMngPtr mng_ptr);
  void init();
  void process();
  void processOneSpectrum(DeconvMsPtr ms_ptr, PrsmWriter &writer);
 private:
    PoissonMngPtr mng_ptr_;
    ProteoformPtrVec proteoforms_;
    PrsmPtrVec prsms_;
    PoissonCompPValuePtr comp_ptr_;
};

typedef std::shared_ptr<PoissonProcessor> PoissonProcessorPtr;

}

#endif
