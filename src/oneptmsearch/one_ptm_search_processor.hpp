#ifndef PROT_ONE_PTM_SEARCH_PROCESSOR_HPP_
#define PROT_ONE_PTM_SEARCH_PROCESSOR_HPP_

#include "htslib/faidx.h"

#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"
#include "oneptmsearch/one_ptm_search_mng.hpp"

namespace prot {

class OnePtmSearchProcessor {
 public:
  OnePtmSearchProcessor(OnePtmSearchMngPtr mng);
  ~OnePtmSearchProcessor();

  void process();

 private:
  OnePtmSearchMngPtr mng_ptr_;
  CompShiftLowMemPtr comp_shift_ptr_;
  // fasta index
  faidx_t *fai_;

  PrsmWriterPtr complete_writer_ptr_;
  PrsmWriterPtr prefix_writer_ptr_;
  PrsmWriterPtr suffix_writer_ptr_;
  PrsmWriterPtr internal_writer_ptr_;

  SimplePrsmReaderPtr complete_reader_ptr_;
  SimplePrsmReaderPtr prefix_reader_ptr_;
  SimplePrsmReaderPtr suffix_reader_ptr_;
  SimplePrsmReaderPtr internal_reader_ptr_;

  ProteoformPtrVec proteo_ptrs_;
  ProteoformPtrVec2D mod_proteo_2d_ptrs_;
  SimplePrsmPtrVec simple_prsm_ptrs_;

  void initReaders();
  void initWriters();
  void initData();
  void closeReaders();
  void closeWriters();
  void processOneSpectrum(SpectrumSetPtr spectrum_set_ptr, 
                          SimplePrsmPtrVec simple_prsm_ptrs,
                          SemiAlignTypePtr align_ptr);
};

typedef std::shared_ptr<OnePtmSearchProcessor> OnePtmSearchProcessorPtr;

} /* namespace prot */

#endif /* ONE_PTM_SEARCH_PROCESSOR_HPP_ */
