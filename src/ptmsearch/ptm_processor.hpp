/*
 * ptm_processor.hpp
 *
 *  Created on: Dec 20, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PTM_PROCESSOR_HPP_
#define PROT_PTM_PROCESSOR_HPP_

#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"
#include "ptmsearch/ptm_mng.hpp"
//#include "ptmsearch/ptm_searcher.hpp"

namespace prot {

class PtmProcessor {
 public:
  PtmProcessor(PtmMngPtr mng);
  void process();

  PtmMngPtr mng_ptr_;
  ProteoformPtrVec proteoforms_;
  SimplePrsmPtrVec simplePrsms_;
  CompShiftLowMemPtr comp_shift_;
  std::vector<PrsmWriterPtr> complete_writers_;
  std::vector<PrsmWriterPtr> prefix_writers_;
  std::vector<PrsmWriterPtr> suffix_writers_;
  std::vector<PrsmWriterPtr> internal_writers_;
  PrsmWriterPtr all_writer_;

 private:
  void init();
  void prsmFindSeq(SimplePrsmPtrVec simple_prsms,ProteoformPtrVec seqs);
  void chooseCompPrePrsms(PrsmPtrVec &all_prsms,
                   PrsmPtrVec &prsms);
  void chooseSuffIntPrsms(PrsmPtrVec &all_prsms,
                     PrsmPtrVec &prsms);
  void search(SpectrumSetPtr spectrum_set_ptr, 
              SimplePrsmPtrVec matches);

};

typedef std::shared_ptr<PtmProcessor> PtmProcessorPtr;

} /* namespace prot */

#endif /* PTM_PROCESSOR_HPP_ */
