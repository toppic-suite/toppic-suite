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

#ifndef TOPPIC_TOPINDEX_PROCESS_HPP
#define TOPPIC_TOPINDEX_PROCESS_HPP

#include <string>
#include <map>
#include <vector>


#include "common/base/mod.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "seq/proteoform.hpp"
#include "seq/db_block.hpp"

#include "filter/zeroptm/mass_zero_ptm_filter.hpp"
#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/diag/diag_filter_mng.hpp"

namespace toppic {
    void TopIndexProcess(std::map<std::string, std::string> & arguments);
    /*
    template <class T>
    T process(T ptr, int thread_num){
        std::string db_file_name = ptr->prsm_para_ptr_->getSearchDbFileName();
        DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

        int block_num = db_block_ptr_vec.size();

        SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);
        
        for (int i = 0; i < block_num; i++) {
            while (pool_ptr->getQueueSize() >= thread_num * 2) {
            boost::this_thread::sleep(boost::posix_time::milliseconds(100));
            }
            pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), ptr));
    }
        pool_ptr->ShutDown();
    };

    template <class T>
    T createIndexFiles(std::string block_str, std::string folderName, T filter_ptr){
        filter_ptr->term_index_ptr_ ->setfileName("term_index" + block_str);
        filter_ptr->diag_index_ptr_->setfileName("diag_index" + block_str);
        filter_ptr->rev_term_index_ptr_->setfileName("rev_term_index" + block_str);
        filter_ptr->rev_diag_index_ptr_->setfileName("rev_diag_index" + block_str);

        filter_ptr->term_index_ptr_->setDirName(folderName);
        filter_ptr->diag_index_ptr_->setDirName(folderName);
        filter_ptr->rev_term_index_ptr_->setDirName(folderName);
        filter_ptr->rev_diag_index_ptr_->setDirName(folderName);
        
        filter_ptr->term_index_ptr_->serializeMassMatch();
        filter_ptr->diag_index_ptr_->serializeMassMatch();
        filter_ptr->rev_term_index_ptr_->serializeMassMatch();
        filter_ptr->rev_diag_index_ptr_->serializeMassMatch();
    };
    */
    void process(PrsmParaPtr prsm_para_ptr, int thread_num);
    //void process(OnePtmFilterMngPtr one_ptm_filter_mng_ptr, int thread_num);
    //void process(DiagFilterMngPtr diag_filter_mng_ptr, int thread_num);
    
    //void createIndexFiles(ProteoformPtrVec &raw_forms, int block_idx, ZeroPtmFilterMngPtr mng_ptr);


}  // namespace toppic

#endif
