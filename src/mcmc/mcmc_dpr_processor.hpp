//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_MCMC_DPR_PROCESSOR_HPP_
#define TOPPIC_MCMC_DPR_PROCESSOR_HPP_


#include <random>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include "seq/proteoform.hpp"
#include "common/base/activation.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "spec/deconv_ms.hpp"

#include "prsm/prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"

#include "tdgf/count_test_num.hpp"
#include "tdgf/comp_pvalue_lookup_table.hpp"

#include "mcmc/mcmc_mng.hpp"

namespace toppic {

class DprProcessor {
 public:
  explicit DprProcessor(MCMCMngPtr mng_ptr):
      mng_ptr_(mng_ptr), mt_(42) {
        init();
      }

  void process();

 private:
  void init();

  std::vector<std::vector<double> > compPtmComb();

  void processOnePrsm(PrsmPtr prsm_ptr, SpectrumSetPtr spec_set_ptr, PrsmXmlWriterPtr prsm_writer);

  MCMCMngPtr mng_ptr_;

  CountTestNumPtr test_num_ptr_;

  std::mt19937 mt_;

  SpParaPtr sp_para_ptr_;

  PtmPtrVec ptm_vec_;

  std::vector<std::vector<double> > ptm_mass_vec2d_;

  std::map<PtmPtr, std::vector<ResiduePtr> > ptm_residue_map_;

  std::map<int, std::vector<std::string> > mass_table_;

  PrsmXmlWriterPtrVec writer_ptr_vec_;

  std::shared_ptr<SimpleThreadPool> pool_ptr_;
};

typedef std::shared_ptr<DprProcessor> DprProcessorPtr;

}  // namespace toppic

#endif
