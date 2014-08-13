#ifndef PROT_PTM_PROCESSOR_HPP_
#define PROT_PTM_PROCESSOR_HPP_

#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"

namespace prot {

class PtmProcessor {
 public:
  PtmProcessor(PtmMngPtr mng);
  void process();

 private:
  PtmMngPtr mng_ptr_;
  ProteoformPtrVec proteoforms_;
  SimplePrsmPtrVec simple_prsms_;
  CompShiftLowMemPtr comp_shift_;
  std::vector<PrsmWriterPtr> complete_writers_;
  std::vector<PrsmWriterPtr> prefix_writers_;
  std::vector<PrsmWriterPtr> suffix_writers_;
  std::vector<PrsmWriterPtr> internal_writers_;
  PrsmWriterPtr all_writer_;

  void initWriters();
  void initData();
  void closeWriters();
  void processOneSpectrum(SpectrumSetPtr spectrum_set_ptr, 
                          SimplePrsmPtrVec matches);
};

typedef std::shared_ptr<PtmProcessor> PtmProcessorPtr;

} /* namespace prot */

#endif /* PTM_PROCESSOR_HPP_ */
