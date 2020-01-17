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

#include "console/topindex_argument.hpp"
#include "console/topindex_process.hpp"

#include "common/base/base_data.hpp"

#include "common/thread/simple_thread_pool.hpp"

#include "seq/fasta_util.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"
#include "prsm/prsm_para.hpp"
#include "seq/proteoform_util.hpp"

#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/diag/diag_filter_mng.hpp"
#include "filter/zeroptm/mass_zero_ptm_filter.hpp"

#include "filter/massmatch/mass_match.hpp"
#include "filter/massmatch/mass_match_factory.hpp"

namespace toppic {
    //create the four pointers and save the data to the output directory
    //total would be 12 files (4 for each type, zero ptm, non ptm, ....)

void TopIndexProcess(std::map<std::string, std::string> &arguments){
    Argument::outputArguments(std::cout, arguments);
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    //so this pointer has all parameter information like thread, database file name. 

    base_data::init();

    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name = ori_db_file_name + "_target";

    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    int thread_num = std::stoi(arguments["threadNumber"]);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }

    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
    
    ZeroPtmFilterMngPtr zero_filter_mng_ptr
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, thread_num, "toppic_zero_filter");

    //OnePtmFilterMngPtr one_ptm_filter_mng_ptr
        //= std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "toppic_one_filter", thread_num);

   // DiagFilterMngPtr diag_filter_mng_ptr
       // = std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num, thread_num, "toppic_multi_filter");

   // process(zero_filter_mng_ptr, thread_num);
    
    //process(one_ptm_filter_mng_ptr, thread_num);
    //process(diag_filter_mng_ptr, thread_num);


}
/*
void createIndexFiles(ProteoformPtrVec &raw_forms, int block_idx, ZeroPtmFilterMngPtr mng_ptr){
    std::string block_str = str_util::toString(block_idx);
    MassZeroPtmFilterPtr filter_ptr = std::make_shared<MassZeroPtmFilter>(raw_forms, mng_ptr, block_str);
    //the pointers containing data are stored inside filter_ptr

    std::string folderName = mng_ptr->prsm_para_ptr_->getOriDbName();
    folderName = folderName + "_index";

    term_index_ptr ->setfileName("term_index" + block_str);
    diag_index_ptr->setfileName("diag_index" + block_str);
    rev_term_index_ptr->setfileName("rev_term_index" + block_str);
    rev_diag_index_ptr->setfileName("rev_diag_index" + block_str);

    term_index_ptr->setDirName(folderName);
    diag_index_ptr->setDirName(folderName);
    rev_term_index_ptr->setDirName(folderName);
    rev_diag_index_ptr->setDirName(folderName);
    
    term_index_ptr->serializeMassMatch();
    diag_index_ptr->serializeMassMatch();
    rev_term_index_ptr->serializeMassMatch();
    rev_diag_index_ptr->serializeMassMatch();
}

std::function<void()> geneTask(int block_idx, ZeroPtmFilterMngPtr mng_ptr){
    return[block_idx, mng_ptr](){

    std::string db_block_file_name = mng_ptr->prsm_para_ptr_->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);    

    ProteoformPtrVec raw_forms 
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name, mng_ptr->prsm_para_ptr_->getFixModPtrVec());

    createIndexFiles(raw_forms, block_idx, mng_ptr);
    };
}

void process(ZeroPtmFilterMngPtr zero_ptr, int thread_num){

    std::string db_file_name = zero_ptr->prsm_para_ptr_->getSearchDbFileName();
    DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

    int block_num = db_block_ptr_vec.size();

    SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);
    
    for (int i = 0; i < block_num; i++) {
        while (pool_ptr->getQueueSize() >= thread_num * 2) {
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
        pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), zero_ptr));
  }
  pool_ptr->ShutDown();
}*/
}