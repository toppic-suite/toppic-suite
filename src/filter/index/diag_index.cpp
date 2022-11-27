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

#include <boost/thread/mutex.hpp>

#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "seq/db_block.hpp"
#include "seq/proteoform_factory.hpp"

#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/index/diag_index.hpp"

namespace toppic{

namespace diag_index {

// serialization mutex.
boost::mutex serial_mutex;

void geneDiagIndexFile(const ProteoformPtrVec &proteo_ptrs,
                       DiagFilterMngPtr mng_ptr, std::string block_str) {

  MassMatchPtr index_ptr 
      = mass_match_factory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                   mng_ptr->max_proteoform_mass_,
                                                   mng_ptr->filter_scale_);
                                                        
  std::string parameters = mng_ptr->getIndexFilePara();
  std::string dir_name = mng_ptr->prsm_para_ptr_->getOriDbName() + "_idx";
  std::string file_name = mng_ptr->multi_ptm_file_vec_[0] + parameters + block_str;
  {
    boost::unique_lock<boost::mutex> lock(serial_mutex);
    index_ptr->serializeMassMatch(file_name, dir_name);
  }
}

std::function<void()> geneIndexTask(int block_idx,
                                    DiagFilterMngPtr mng_ptr){
  return [block_idx, mng_ptr] () {

    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder()
        + "_" + str_util::toString(block_idx);
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    std::string block_str = str_util::toString(block_idx);
    geneDiagIndexFile(raw_forms, mng_ptr, block_str);
  };                                      
}

void process(DiagFilterMngPtr mng_ptr) {
  std::string db_file_name = mng_ptr->prsm_para_ptr_->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr->thread_num_);
  int block_num = db_block_ptr_vec.size();
  std::cout << "Generating multiple shift index files --- started" << std::endl;
  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }

    std::cout << "Multiple shift index files - processing " << (i+1) 
        << " of " << block_num << " files." << std::endl;
    pool_ptr->Enqueue(geneIndexTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr));
  }
  pool_ptr->ShutDown();
  std::cout << "Generating multiple shift index files --- finished" << std::endl;
}

}

}
