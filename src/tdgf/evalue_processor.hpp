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


#ifndef PROT_EVALUE_PROCESSOR_HPP_
#define PROT_EVALUE_PROCESSOR_HPP_

#include "htslib/faidx.h"

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_pvalue_array.hpp"
#include "tdgf/comp_pvalue_lookup_table.hpp"

namespace prot {

class EValueProcessor {
 public:
  EValueProcessor(TdgfMngPtr mng_ptr): mng_ptr_(mng_ptr) {};
  ~EValueProcessor();

  void init();
  void process(bool is_separate);
  void processOneSpectrum(SpectrumSetPtr spec_set_ptr,
                          PrsmPtrVec &sele_prsm_ptrs,
                          double ppo, bool is_separate,
                          PrsmXmlWriter &writer);

  void compEvalues(SpectrumSetPtr spectrum_set_ptr, 
                   PrsmPtrVec &sele_prsm_ptrs,
                   double ppo, bool is_separate);

 private:
  TdgfMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;
  CompPValueArrayPtr comp_pvalue_ptr_;
  CompPValueLookupTablePtr comp_pvalue_table_ptr_;
  faidx_t *fai_;
  bool checkPrsms(const PrsmPtrVec &prsm_ptrs);

};

typedef std::shared_ptr<EValueProcessor> EValueProcessorPtr;

}

#endif
