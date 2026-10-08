// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "filter/zeroptm/zero_ptm_filter_processor.hpp"

#include <iomanip>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "filter/massmatch/filter_spec.hpp"
#include "filter/zeroptm/zero_ptm_filter.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/simple_prsm_mem_merge.hpp"
#include "prsm/simple_prsm_str_merge.hpp"
#include "prsm/simple_prsm_xml_writer_set.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace zero_ptm_filter_processor {

// Reads the msalign file once and converts every spectrum set into the block
// independent form the filter scores, so that the blocks share the conversion
// instead of repeating it.
FilterSpecPtrVec readSpectra(const ZeroPtmFilterMngPtr& mng_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReaderPtr reader_ptr = std::make_shared<MsAlignReader>(
      prsm_para_ptr->getSpectrumFileName(), group_spec_num,
      sp_para_ptr->getActivationPtr(), sp_para_ptr->getEnvCnnCutoff(),
        sp_para_ptr->getFragDedupPpo());
  std::vector<double> prec_error_vec =
      sp_para_ptr->getZeroShiftSearchPrecErrorVec();

  FilterSpecPtrVec spectra;
  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  while (deconv_ms_ptr_vec.size() != 0) {
    if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->containsPrec()) {
      // allow one dalton error
      SpectrumSetPtrVec spec_set_vec =
          spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError(
              deconv_ms_ptr_vec, sp_para_ptr, prec_error_vec);
      for (size_t k = 0; k < spec_set_vec.size(); k++) {
        if (spec_set_vec[k]->isValid()) {
          spectra.push_back(ZeroPtmFilter::geneFilterSpec(
              spec_set_vec[k]->getMsThreePtrVec(), mng_ptr));
        }
      }
    }
    deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  }
  return spectra;
}

void scoreBlock(const ZeroPtmFilterPtr& filter_ptr, int block_idx,
                const ZeroPtmFilterMngPtr& mng_ptr,
                const FilterSpecPtrVec& spectra, SimplePrsmMemMerge* combiner) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  int block_num = static_cast<int>(mng_ptr->cnts_.size());
  int spec_num = mng_ptr->n_spec_block_ / block_num;
  size_t report_step = 100;
  for (size_t i = 0; i < spectra.size(); i++) {
    filter_ptr->computeBestMatch(spectra[i]);
    int add_idx = static_cast<int>(i);
    combiner->add(0, block_idx, add_idx, filter_ptr->getCompMatchPtrs());
    combiner->add(1, block_idx, add_idx, filter_ptr->getPrefMatchPtrs());
    combiner->add(2, block_idx, add_idx, filter_ptr->getSuffMatchPtrs());
    combiner->add(3, block_idx, add_idx, filter_ptr->getInternalMatchPtrs());

    if ((i + 1) % report_step == 0 || i + 1 == spectra.size()) {
      // the count covers the spectra without precursor as well
      mng_ptr->cnts_[block_idx] = i + 1 == spectra.size()
                                      ? spec_num
                                      : static_cast<int>((i + 1) *
                                                         group_spec_num);
      int cnt_sum = 0;
      for (size_t k = 0; k < mng_ptr->cnts_.size(); k++) {
        cnt_sum = cnt_sum + mng_ptr->cnts_[k];
      }
      double perc = cnt_sum * 100.0 / mng_ptr->n_spec_block_;
      std::stringstream msg;
      msg << std::flush << "Zero unexpected shift filtering - processing "
          << std::setprecision(3) << perc << "%.    \r";
      mng_ptr->mutex_.lock();
      std::cout << msg.str();
      mng_ptr->mutex_.unlock();
    }
  }
}

std::function<void()> geneTask(int block_idx,
                               const ZeroPtmFilterMngPtr& mng_ptr,
                               const FilterSpecPtrVec* spectra,
                               SimplePrsmMemMerge* combiner) {
  return [block_idx, mng_ptr, spectra, combiner]() {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string db_block_file_name =
        prsm_para_ptr->getSearchDbFileNameWithFolder() + "_" +
        std::to_string(block_idx);

    ProteoformPtrVec raw_forms =
        proteoform_factory::readFastaToProteoformPtrVec(
            db_block_file_name, prsm_para_ptr->getFixModPtrVec());
    ZeroPtmFilterPtr filter_ptr = std::make_shared<ZeroPtmFilter>(
        raw_forms, mng_ptr, std::to_string(block_idx));
    scoreBlock(filter_ptr, block_idx, mng_ptr, *spectra, combiner);
  };
}

void process(const ZeroPtmFilterMngPtr& mng_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  int spec_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());

  mng_ptr->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  // the candidates of the blocks are combined in memory instead of in one
  // file per block and proteoform type
  std::vector<int> top_nums = {static_cast<int>(mng_ptr->comp_num_),
                               static_cast<int>(mng_ptr->pref_suff_num_),
                               static_cast<int>(mng_ptr->pref_suff_num_),
                               static_cast<int>(mng_ptr->inte_num_)};
  SimplePrsmMemMerge combiner(
      prsm_para_ptr->getSpectrumFileName(), mng_ptr->output_file_ext_,
      top_nums, msalign_util::readSpIds(prsm_para_ptr->getSpectrumFileName()));
  int block_num = db_block_ptr_vec.size();
  mng_ptr->cnts_.resize(block_num, 0);
  // the spectra are read and converted once and scored against every block
  FilterSpecPtrVec spectra = readSpectra(mng_ptr);
  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(mng_ptr->thread_num_);

  LOG_DEBUG("thread num " << mng_ptr->thread_num_);
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    pool_ptr->enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr,
                               &spectra, &combiner));
  }
  pool_ptr->shutDown();
  std::cout << std::endl;
  std::cout << "Zero unexpected shift filtering - combining blocks started."
            << std::endl;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_pref = mng_ptr->output_file_ext_;
  combiner.write();
  std::cout << "Zero unexpected shift filtering - combining blocks finished."
            << std::endl;
}

}  // namespace zero_ptm_filter_processor

}  // namespace toppic
