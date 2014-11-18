#ifndef PROT_PTM_PROCESSOR_HPP_
#define PROT_PTM_PROCESSOR_HPP_

#include "htslib/faidx.h"

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
  ~PtmProcessor();

  void process();

 private:
  PtmMngPtr mng_ptr_;
  CompShiftLowMemPtr comp_shift_ptr_;
  // fasta index
  faidx_t *fai_;

  std::vector<PrsmWriterPtr> complete_writer_ptrs_;
  std::vector<PrsmWriterPtr> prefix_writer_ptrs_;
  std::vector<PrsmWriterPtr> suffix_writer_ptrs_;
  std::vector<PrsmWriterPtr> internal_writer_ptrs_;
  PrsmWriterPtr all_writer_ptr_;

  ProteoformPtrVec proteo_ptrs_;
  ProteoformPtrVec2D mod_proteo_2d_ptrs_;
  SimplePrsmPtrVec simple_prsm_ptrs_;

  void initWriters();
  void initData();
  void closeWriters();
  void processOneSpectrum(SpectrumSetPtr spectrum_set_ptr, 
                          SimplePrsmPtrVec simple_prsm_ptrs);
};

typedef std::shared_ptr<PtmProcessor> PtmProcessorPtr;

} /* namespace prot */

#endif /* PTM_PROCESSOR_HPP_ */
