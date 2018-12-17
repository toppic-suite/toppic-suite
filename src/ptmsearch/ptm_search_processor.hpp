//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef PROT_PTM_SEARCH_PROCESSOR_HPP_
#define PROT_PTM_SEARCH_PROCESSOR_HPP_

#include "base/fasta_index_reader.hpp"
#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"

namespace toppic {

template <int N>
class PrsmXmlWriterSet {
 public:
  PrsmXmlWriterSet(const std::string & output_file_name);
  std::vector<PrsmXmlWriterPtr> complete_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> prefix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> suffix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> internal_writer_ptrs_;
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

} /* namespace toppic */

#endif /* PTM_PROCESSOR_HPP_ */
