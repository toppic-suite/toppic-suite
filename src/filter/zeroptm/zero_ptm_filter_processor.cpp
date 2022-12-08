//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <iomanip>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

#include "ms/spec/msalign_util.hpp"
#include "ms/factory/spectrum_set_factory.hpp"

#include "prsm/simple_prsm_xml_writer_set.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

#include "filter/zeroptm/zero_ptm_filter.hpp"
#include "filter/zeroptm/zero_ptm_filter_processor.hpp"

namespace toppic {

namespace zero_ptm_filter_processor {

inline void filterBlock(const ProteoformPtrVec & raw_forms,
                        int block_idx, ZeroPtmFilterMngPtr mng_ptr) { 
  std::string block_str = str_util::toString(block_idx);
  ZeroPtmFilterPtr filter_ptr = std::make_shared<ZeroPtmFilter>(raw_forms, mng_ptr, block_str);

  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  SimpleMsAlignReaderPtr reader_ptr 
      = std::make_shared<SimpleMsAlignReader>(prsm_para_ptr->getSpectrumFileName(), 
                                              group_spec_num,
                                              sp_para_ptr->getActivationPtr());
  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr->output_file_ext_ + "_" + block_str;
  //writer
  SimplePrsmXmlWriterSet writers(output_file_name);
  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  std::vector<double> prec_error_vec = sp_para_ptr->getZeroShiftSearchPrecErrorVec(); 
  while (deconv_ms_ptr_vec.size() != 0) {
    // allow one dalton error
    SpectrumSetPtrVec spec_set_vec 
        = spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError(deconv_ms_ptr_vec, 
                                                                   sp_para_ptr,
                                                                   prec_error_vec);
    for (size_t k = 0; k < spec_set_vec.size(); k++) {
      LOG_DEBUG("spec set ptr valid " << spec_set_vec[k]->isValid());
      if (spec_set_vec[k]->isValid()) {
        ExtendMsPtrVec ms_ptr_vec = spec_set_vec[k]->getMsThreePtrVec();
        filter_ptr->computeBestMatch(ms_ptr_vec);
        writers.getCompleteWriterPtr()->write(filter_ptr->getCompMatchPtrs());
        writers.getPrefixWriterPtr()->write(filter_ptr->getPrefMatchPtrs());
        writers.getSuffixWriterPtr()->write(filter_ptr->getSuffMatchPtrs());
        writers.getInternalWriterPtr()->write(filter_ptr->getInternalMatchPtrs());
      }
    }

    mng_ptr->cnts_[block_idx] = mng_ptr->cnts_[block_idx] + group_spec_num;
    int cnt_sum = 0;
    for (size_t i = 0; i < mng_ptr->cnts_.size(); i++) {
      cnt_sum = cnt_sum + mng_ptr->cnts_[i];
    }
    double perc = cnt_sum * 100.0 / mng_ptr->n_spec_block_;
    std::stringstream msg;
    msg << std::flush << "Zero unexpected shift filtering - processing " << std::setprecision(3) <<  perc << "%.    \r";
    mng_ptr->mutex_.lock();
    std::cout << msg.str();
    mng_ptr->mutex_.unlock();
    deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  }
  writers.close();
}

std::function<void()> geneTask(int block_idx, 
                               ZeroPtmFilterMngPtr mng_ptr) {
  return[block_idx, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder()
      + "_" + str_util::toString(block_idx);

    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());

    filterBlock(raw_forms, block_idx, mng_ptr);
  };
}

void process(ZeroPtmFilterMngPtr mng_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  int spec_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());

  mng_ptr->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr->thread_num_);
  int block_num = db_block_ptr_vec.size();
  mng_ptr->cnts_.resize(block_num, 0);

  LOG_DEBUG("thread num " << mng_ptr->thread_num_);
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() ==0) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr));
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;
  std::cout << "Zero unexpected shift filtering - combining blocks started." << std::endl;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_pref = mng_ptr->output_file_ext_;
  SimplePrsmStrMerge::mergeBlockResults(sp_file_name, input_pref, block_num,  
                                        mng_ptr->comp_num_, mng_ptr->pref_suff_num_, mng_ptr->inte_num_ );
  std::cout << "Zero unexpected shift filtering - combining blocks finished." << std::endl;
}

}

}  // namespace toppic
