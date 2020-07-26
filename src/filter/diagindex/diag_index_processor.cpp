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

#include <iostream>

#include "common/util/file_util.hpp"
#include "common/base/mod_util.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "seq/proteoform_factory.hpp"

#include "filter/mng/topindex_file_name.hpp"

#include "filter/diagindex/diag_index_processor.hpp"
#include "filter/diagindex/diag_index_file.hpp"

namespace toppic{

std::function<void()> geneIndexTask(int block_idx,
                                    DiagFilterMngPtr mng_ptr){
return [block_idx, mng_ptr] () {

  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
      + "_" + str_util::toString(block_idx);
  ProteoformPtrVec raw_forms
      = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                        prsm_para_ptr->getFixModPtrVec());
  std::string block_str = str_util::toString(block_idx);
  diag_index_file::geneDiagIndexFile(raw_forms, mng_ptr, block_str);
  };                                      
}

void DiagIndexProcessor::process() {
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  std::cout << "Generating multiple PTM index files --- started" << std::endl;
  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    
    std::cout << "Multiple PTM index files - processing " << (i+1) << " of " << block_num << " files." << std::endl;
    pool_ptr->Enqueue(geneIndexTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr_));
  }
  pool_ptr->ShutDown();
  std::cout << "Generating multiple PTM index files --- finished" << std::endl;
}



}
