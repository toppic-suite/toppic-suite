#ifndef PROT_PTM_SEARCH_PROCESSOR_HPP_
#define PROT_PTM_SEARCH_PROCESSOR_HPP_

#include "base/fasta_index_reader.hpp"
#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"

namespace prot {

class PtmSearchProcessor {
 public:
  PtmSearchProcessor(PtmSearchMngPtr mng_ptr);

  void process();

 private:
  PtmSearchMngPtr mng_ptr_;
  CompShiftLowMemPtr comp_shift_ptr_;

  std::vector<PrsmXmlWriterPtr> complete_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> prefix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> suffix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> internal_writer_ptrs_;
  PrsmXmlWriterPtr all_writer_ptr_;

  ProteoformPtrVec proteo_ptrs_;
  ProteoformPtrVec2D mod_proteo_2d_ptrs_;
  SimplePrsmPtrVec simple_prsm_ptrs_;

  void initWriters();
  void closeWriters();
  void processOneSpectrum(SpectrumSetPtr spectrum_set_ptr, 
                          SimplePrsmPtrVec simple_prsm_ptrs);
};

typedef std::shared_ptr<PtmSearchProcessor> PtmSearchProcessorPtr;

} /* namespace prot */

#endif /* PTM_PROCESSOR_HPP_ */
