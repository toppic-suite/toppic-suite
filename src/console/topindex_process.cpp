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

#include "common/base/mod.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/base_data.hpp"

#include "common/thread/simple_thread_pool.hpp"

#include "seq/fasta_util.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "filter/diag/diag_filter_mng.hpp"
#include "filter/diag/diag_filter_processor.hpp"

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

    std::vector<std::string> index_file_names;

    //prsm_para_ptr->fixed_mod_list = mod_util::geneFixedModList(arguments["fixedMod"]);

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }

    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
    
    ZeroPtmFilterMngPtr zero_filter_mng_ptr
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, thread_num, "toppic_zero_filter");

    //ZeroPtmFilterProcessorPtr zero_filter_processor
       // = std::make_shared<ZeroPtmFilterProcessor>(zero_filter_mng_ptr);

    OnePtmFilterMngPtr one_ptm_filter_mng_ptr
        = std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "toppic_one_filter", thread_num);
    //OnePtmFilterProcessorPtr one_filter_processor
        //= std::make_shared<OnePtmFilterProcessor>(one_ptm_filter_mng_ptr);

    DiagFilterMngPtr diag_filter_mng_ptr
        = std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num, thread_num, "toppic_multi_filter");
    //DiagFilterProcessorPtr diag_filter_processor
        //= std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);

    index_file_names.push_back("toppic_zero_ptm_complete");
    index_file_names.push_back("toppic_zero_ptm_prefix");
    index_file_names.push_back("toppic_zero_ptm_suffix");
    index_file_names.push_back("toppic_zero_ptm_internal");  

    index_file_names.push_back("toppic_one_ptm_complete");
    index_file_names.push_back("toppic_one_ptm_prefix");
    index_file_names.push_back("toppic_one_ptm_suffix");
    index_file_names.push_back("toppic_one_ptm_internal");

    index_file_names.push_back("toppic_multi_ptm_complete_2");
    index_file_names.push_back("toppic_multi_ptm_prefix_2");
    index_file_names.push_back("toppic_multi_ptm_suffix_2");
    index_file_names.push_back("toppic_multi_ptm_internal_2");

    zero_filter_mng_ptr = prsm_para_ptr;
    zero_filter_mng_ptr->thread_num_ = thread_num;

    one_ptm_filter_mng_ptr = prsm_para_ptr;
    one_ptm_filter_mng_ptr->thread_num_ = thread_num;

    diag_filter_mng_ptr = prsm_para_ptr;
    diag_filter_mng_ptr->thread_num_ = thread_num;

    //process(prsm_para_ptr, zero_filter_processor, one_filter_processor, diag_filter_processor, thread_num);
    //maybe should I process each filter processors one by one..

}

void createIndexFiles(ProteoformPtrVec raw_forms, int block_idx, ZeroPtmFilterMngPtr mng_ptr){
    std::string block_str = str_util::toString(block_idx);
    MassZeroPtmFilterPtr filter_ptr = std::make_shared<MassZeroPtmFilter>(raw_forms, mng_ptr, block_str);
    //create the index files here, and restore the mass zero ptm filter file to previous version (without serialization)
}

std::function<void()> geneTask(int block_idx, ZeroPtmFilterMngPtr mng_ptr){
    return[block_idx, mng_ptr](){

    std::string db_block_file_name = mng_ptr->prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);    

    ProteoformPtrVec raw_forms 
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name, prsm_para_ptr->getFixModPtrVec());

    createIndexFiles(raw_forms, block_idx, mng_ptr);
    };
}

void process(ZeroPtmFilterProcessorPtr zero_ptr){

    std::string db_file_name = zero_ptr->prsm_para_ptr->getSearchDbFileName();
    DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

    int block_num = db_block_ptr_vec.size();

    SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(zero_ptr->thread_num);
    
    for (int i = 0; i < block_num; i++) {
        while (pool_ptr->getQueueSize() >= zero_ptr->thread_num * 2) {
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
        pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), zero_ptr));
  }
  pool_ptr->ShutDown();


}
}