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
  SimplePrSMPtrVec simplePrsms_;
  CompShiftLowMemPtr comp_shift_;
  std::vector<PrSMWriterPtr> complete_writers_;
  std::vector<PrSMWriterPtr> prefix_writers_;
  std::vector<PrSMWriterPtr> suffix_writers_;
  std::vector<PrSMWriterPtr> internal_writers_;
  PrSMWriterPtr all_writer_;

 private:
  void init();
  void prsmFindSeq(SimplePrSMPtrVec simple_prsms,ProteoformPtrVec seqs);
  void chooseCompPrePrsms(PrSMPtrVec &all_prsms,
                   PrSMPtrVec &prsms);
  void chooseSuffIntPrsms(PrSMPtrVec &all_prsms,
                     PrSMPtrVec &prsms);
  void search(SpectrumSetPtr spectrum_set_ptr, 
              SimplePrSMPtrVec matches
              );

};

typedef std::shared_ptr<PtmProcessor> PtmProcessorPtr;

} /* namespace prot */

#endif /* PTM_PROCESSOR_HPP_ */
