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
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/spectrum_set.hpp"

#include "prsm/simple_prsm_xml_writer_set.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "filter/zeroptm/mass_zero_ptm_filter.hpp"

#include "filter/zeroptm/mass_zero_ptm_index_file.hpp"
#include "filter/zeroptm/zero_ptm_count_mng.hpp"

#include "console/topindex_file_name.hpp"

namespace toppic {

inline void filterBlock(const ProteoformPtrVec & raw_forms,
                        int block_idx, ZeroPtmFilterMngPtr mng_ptr) { 
  std::string block_str = str_util::toString(block_idx);
  int group_spec_num = mng_ptr->prsm_para_ptr_->getGroupSpecNum();
  MassZeroPtmFilterPtr filter_ptr = std::make_shared<MassZeroPtmFilter>(raw_forms, mng_ptr, block_str);
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       prsm_para_ptr->getSpParaPtr()->getActivationPtr(),
                       prsm_para_ptr->getSpParaPtr()->getSkipList());
  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr->output_file_ext_ + "_" + block_str;
  //writer
  SimplePrsmXmlWriterSet writers(output_file_name);
  std::vector<SpectrumSetPtr> spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
  while (spec_set_vec[0] != nullptr) {

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
    mng_ptr->cnts_[block_idx] = mng_ptr->cnts_[block_idx] + 1;
    int cnt_sum = 0;
    for (size_t i = 0; i < mng_ptr->cnts_.size(); i++) {
      cnt_sum = cnt_sum + mng_ptr->cnts_[i];
    }
    double perc = cnt_sum * 100.0 / mng_ptr->n_spec_block_;
    std::stringstream msg;
    msg << std::flush << "Non PTM filtering - processing " << std::setprecision(3) <<  perc << "%.    \r";
    mng_ptr->mutex_.lock();
    std::cout << msg.str();
    mng_ptr->mutex_.unlock();
    spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
  }
  reader.close();
  writers.close();
}

std::function<void()> geneTask(int block_idx, 
                               ZeroPtmFilterMngPtr mng_ptr) {
  return[block_idx, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    //std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);

    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());

    filterBlock(raw_forms, block_idx, mng_ptr);
  };
}

void ZeroPtmFilterProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  int spec_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  // n_spec_block = spec_num * block_num
  mng_ptr_->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  mng_ptr_->cnts_.resize(block_num, 0);
  //logger::setLogLevel(2);
  LOG_DEBUG("thread num " << mng_ptr_->thread_num_);

  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr_));
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;
  std::cout << "Non PTM filtering - combining blocks started." << std::endl;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_pref = mng_ptr_->output_file_ext_;
  SimplePrsmStrMerge::mergeBlockResults(sp_file_name, input_pref, block_num,  
                                        mng_ptr_->comp_num_, mng_ptr_->pref_suff_num_, mng_ptr_->inte_num_ );
  std::cout << "Non PTM filtering - combining blocks finished." << std::endl;
}

//below functions are used for generating index files

inline void createIndexFiles(const ProteoformPtrVec & raw_forms,
                          int block_idx, ZeroPtmFilterMngPtr mng_ptr, int block_num, int *current_num) {  
                        
    std::string block_str = str_util::toString(block_idx);
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    //index file names
    //naming = fixed mod + n terminal + activation + error tolerance + decoy + block number
    
    TopIndexFileName TopIndexFile;
    std::string parameters = TopIndexFile.gene_file_name(prsm_para_ptr);

    std::vector<std::string> file_vec{TopIndexFile.zero_ptm_file_vec[0] + parameters + block_str, 
    TopIndexFile.zero_ptm_file_vec[1] + parameters + block_str, 
    TopIndexFile.zero_ptm_file_vec[2] + parameters + block_str, TopIndexFile.zero_ptm_file_vec[3] + parameters + block_str};

    MassZeroPtmIndexPtr filter_ptr = std::make_shared<MassZeroPtmIndex>(raw_forms, mng_ptr, file_vec);

    mng_ptr->mutex_.lock();

    std::cout << "Non PTM index files - processing " << *current_num << " of " << block_num << " files." << std::endl;
    *current_num = *current_num + 1;

    mng_ptr->mutex_.unlock();
  
}

std::function<void()> geneIndexTask(int block_idx, 
                               ZeroPtmFilterMngPtr mng_ptr, int block_num, int *current_num) {
  return[block_idx, mng_ptr, block_num, current_num] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());

    createIndexFiles(raw_forms, block_idx, mng_ptr, block_num, current_num);
  };
}
void ZeroPtmFilterProcessor::index_process(){
  //for generating index files
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  // n_spec_block = spec_num * block_num

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  
  int block_num = db_block_ptr_vec.size();
  int current_num = 1; //show how many files have been processed. n in the message "n of 5 files processed"..

  //logger::setLogLevel(2);
  std::cout << "Generating Non PTM index files --- started" << std::endl;
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneIndexTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr_, block_num, &current_num));
  }
  pool_ptr->ShutDown();
  std::cout << "Generating Non PTM index files --- finished" << std::endl;
  //std::cout << std::endl;
}

}  // namespace toppic
