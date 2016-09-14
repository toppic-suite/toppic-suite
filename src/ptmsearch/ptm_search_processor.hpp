// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
