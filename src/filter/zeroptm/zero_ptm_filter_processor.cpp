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

#include <iomanip>

#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/spectrum_set.hpp"

#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"

#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "filter/zeroptm/mass_zero_ptm_filter.hpp"

namespace toppic {

void ZeroPtmFilterProcessor::combineBlockResults() {
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);
  std::string output_pref = mng_ptr_->output_file_ext_ + "_";

  std::string complete = ProteoformType::COMPLETE->getName();
  std::string prefix = ProteoformType::PREFIX->getName();
  std::string suffix = ProteoformType::SUFFIX->getName();
  std::string internal = ProteoformType::INTERNAL->getName();

  int block_num = db_block_ptr_vec.size();

  SimplePrsmStrCombine comp_combine(sp_file_name, 
                                    output_pref + complete,
                                    block_num, 
                                    output_pref + complete,
                                    mng_ptr_->comp_num_);
  comp_combine.process();

  SimplePrsmStrCombine pref_combine(sp_file_name, 
                                    output_pref + prefix,
                                    block_num, 
                                    output_pref + prefix,
                                    mng_ptr_->pref_suff_num_);
  pref_combine.process();

  SimplePrsmStrCombine suff_combine(sp_file_name, 
                                    output_pref + suffix,
                                    block_num, 
                                    output_pref + suffix,
                                    mng_ptr_->pref_suff_num_);
  suff_combine.process();

  SimplePrsmStrCombine internal_combine(sp_file_name, 
                                        output_pref + internal,
                                        block_num, 
                                        output_pref + internal,
                                        mng_ptr_->inte_num_);
  internal_combine.process();

  // remove tempory files
  std::string end_str = "_";
  file_util::cleanTempFiles(sp_file_name, output_pref + complete + end_str);
  file_util::cleanTempFiles(sp_file_name, output_pref + prefix + end_str);
  file_util::cleanTempFiles(sp_file_name, output_pref + suffix + end_str);
  file_util::cleanTempFiles(sp_file_name, output_pref + internal + end_str);
}

inline void filterBlock(const ProteoformPtrVec & raw_forms,
                        const std::string & block_str,
                        ZeroPtmFilterMngPtr mng_ptr) {
  int group_spec_num = mng_ptr->prsm_para_ptr_->getGroupSpecNum();
  MassZeroPtmFilterPtr filter_ptr = std::make_shared<MassZeroPtmFilter>(raw_forms, mng_ptr);
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       prsm_para_ptr->getSpParaPtr()->getActivationPtr(),
                       prsm_para_ptr->getSpParaPtr()->getSkipList());
  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr->output_file_ext_;
  std::string complete = ProteoformType::COMPLETE->getName();
  std::string prefix = ProteoformType::PREFIX->getName();
  std::string suffix = ProteoformType::SUFFIX->getName();
  std::string internal = ProteoformType::INTERNAL->getName();
  SimplePrsmXmlWriter comp_writer(output_file_name + "_" + complete + "_" + block_str);
  SimplePrsmXmlWriter pref_writer(output_file_name + "_" + prefix + "_" + block_str);
  SimplePrsmXmlWriter suff_writer(output_file_name + "_" + suffix + "_" + block_str);
  SimplePrsmXmlWriter internal_writer(output_file_name + "_" + internal + "_" + block_str);

  std::vector<SpectrumSetPtr> spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
  while (spec_set_vec[0] != nullptr) {
    for (size_t k = 0; k < spec_set_vec.size(); k++) {
      LOG_DEBUG("spec set ptr valid " << spec_set_vec[k]->isValid());
      if (spec_set_vec[k]->isValid()) {
        ExtendMsPtrVec ms_ptr_vec = spec_set_vec[k]->getMsThreePtrVec();
        filter_ptr->computeBestMatch(ms_ptr_vec);
        comp_writer.write(filter_ptr->getCompMatchPtrs());
        pref_writer.write(filter_ptr->getPrefMatchPtrs());
        suff_writer.write(filter_ptr->getSuffMatchPtrs());
        internal_writer.write(filter_ptr->getInternalMatchPtrs());
      }
    }
    mng_ptr->cnt_++;
    double perc = mng_ptr->cnt_ * 100.0 / mng_ptr->n_spec_block_;
    std::cout << std::flush << "Non PTM filtering - processing " << std::setprecision(3) <<  perc << "%.\r";
    spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
  }
  reader.close();
  comp_writer.close();
  pref_writer.close();
  suff_writer.close();
  internal_writer.close();
}

std::function<void()> geneTask(DbBlockPtr block_ptr, 
                               ZeroPtmFilterMngPtr mng_ptr) {
  return[block_ptr, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_ptr->getBlockIdx());
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    std::string block_str = str_util::toString(block_ptr->getBlockIdx());
    filterBlock(raw_forms, block_str, mng_ptr);
  };
}

void ZeroPtmFilterProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::cout << "Non PTM filtering started " << std::endl;
  // cnt_ is thread_safe 
  mng_ptr_->cnt_ = 0;
  int spec_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  // n_spec_block = spec_num * block_num
  mng_ptr_->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i], mng_ptr_));
  }
  pool_ptr->ShutDown();
  std::cout << "Non PTM filtering finished. " << std::endl;

  std::cout << "Non PTM filtering - combining blocks started." << std::endl;
  combineBlockResults();
  std::cout << "Non PTM filtering - combining blocks finished." << std::endl;
}

}  // namespace toppic
