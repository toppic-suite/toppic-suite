#ifndef PROT_POSSION_PROCESSOR_HPP_
#define PROT_POSSION_PROCESSOR_HPP_

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "possion/possion_mng.hpp"
#include "possion/possion_comp_pvalue.hpp"

namespace prot {

class PossionProcessor {
 public:
  PossionProcessor(PossionMngPtr mng_ptr);
  void init();
  void process();
  void processOneSpectrum(DeconvMsPtr ms_ptr, PrsmWriter &writer);
 private:
    PossionMngPtr mng_ptr_;
    ProteoformPtrVec proteoforms_;
    PrsmPtrVec prsms_;
    PossionCompPValuePtr comp_ptr_;
};

}

#endif
