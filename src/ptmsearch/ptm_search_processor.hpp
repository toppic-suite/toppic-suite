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

template <int N>
class PrsmXmlWriterSet {
 public:
  PrsmXmlWriterSet(const std::string & output_file_name);
  std::vector<PrsmXmlWriterPtr> complete_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> prefix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> suffix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> internal_writer_ptrs_;
  PrsmXmlWriterPtr all_writer_ptr_;
  void close();
};

class PtmSearchProcessor {
 public:
  PtmSearchProcessor(PtmSearchMngPtr mng_ptr);

  void process();

 private:
  PtmSearchMngPtr mng_ptr_;
  CompShiftLowMem comp_shift_;

  ProteoformPtrVec proteo_ptrs_;
  ProteoformPtrVec2D mod_proteo_2d_ptrs_;
  SimplePrsmPtrVec simple_prsm_ptrs_;

};

typedef std::shared_ptr<PtmSearchProcessor> PtmSearchProcessorPtr;

} /* namespace prot */

#endif /* PTM_PROCESSOR_HPP_ */
