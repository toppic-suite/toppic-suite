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
#include "seq/proteoform_util.hpp"

#include "prsm/prsm_para.hpp"

#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/zeroptm/mass_zero_ptm_filter.hpp"

#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/oneptm/mass_one_ptm_filter.hpp"

#include "filter/diag/diag_filter_mng.hpp"
#include "filter/diag/mass_diag_filter.hpp"

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
    
    process(prsm_para_ptr, thread_num);


    //ZeroPtmFilterMngPtr zero_ptm_index_ptr = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr);

   // OnePtmFilterMngPtr one_ptm_filter_index_ptr = std::make_shared<OnePtmFilterMng>(prsm_para_ptr);

    //DiagFilterMngPtr diag_filter_index_ptr = std::make_shared<DiagFilterMng>(prsm_para_ptr);

    //
    //process<OnePtmFilterMngPtr>(one_ptm_filter_index_ptr, thread_num);
    //process<DiagFilterMngPtr>(diag_filter_index_ptr, thread_num);

    //process(one_ptm_filter_mng_ptr, thread_num);
    //process(diag_filter_mng_ptr, thread_num);

}

void createIndexPtr(ProteoformPtrVec &proteo_ptrs, int block_idx, PrsmParaPtr prsm_para_ptr){
    std::string block_str = str_util::toString(block_idx);
    std::string folderName = prsm_para_ptr->getOriDbName();

    ZeroPtmFilterMngPtr zero_ptm_ptr = std::make_shared<ZeroPtmFilterMng>();
    ZeroPtmFilterMngPtr

    folderName = folderName + "_index";

    std::vector<std::vector<double> > shift_2d
      = proteoform_util::getNTermShift2D(proteo_ptrs, prsm_para_ptr->getProtModPtrVec());

    std::vector<std::vector<double> > n_term_acet_2d
      = proteoform_util::getNTermAcet2D(proteo_ptrs, prsm_para_ptr->getProtModPtrVec());

    MassMatchPtr t_ptr = MassMatchFactory::getPrmTermMassMatchPtr(proteo_ptrs, shift_2d,
                                                             zero_ptm_ptr->max_proteoform_mass_,
                                                             zero_ptm_ptr->filter_scale_);

    MassMatchPtr d_ptr = MassMatchFactory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                             zero_ptm_ptr->max_proteoform_mass_,
                                                             zero_ptm_ptr->filter_scale_);

    std::vector<std::vector<double> > rev_shift_2d;
    std::vector<double> shift_1d(1, 0);
    for (size_t i = 0; i < proteo_ptrs.size(); i++) {
        rev_shift_2d.push_back(shift_1d);
    }

    MassMatchPtr rev_t_ptr = MassMatchFactory::getSrmTermMassMatchPtr(proteo_ptrs, rev_shift_2d,
                                                                 n_term_acet_2d,
                                                                 zero_ptm_ptr->max_proteoform_mass_,
                                                                 zero_ptm_ptr->filter_scale_);

    MassMatchPtr rev_d_ptr = MassMatchFactory::getSrmDiagMassMatchPtr(proteo_ptrs, n_term_acet_2d,
                                                                 zero_ptm_ptr->max_proteoform_mass_,
                                                                 zero_ptm_ptr->filter_scale_);                                                        
        
    //t_ptr->serializeMassMatch("term_index", block_str, folderName);
    //d_ptr->serializeMassMatch("diag_index", block_str, folderName);
    //rev_t_ptr->serializeMassMatch("rev_term_index", block_str, folderName);
    //rev_d_ptr->serializeMassMatch("rev_diag_index", block_str, folderName);
}
/*
void createIndexPtr(ProteoformPtrVec &raw_forms, int block_idx, OnePtmFilterMngPtr mng_ptr){
    std::string block_str = str_util::toString(block_idx);
    MassOnePtmFilterPtr filter_ptr = std::make_shared<OnePtmFilter>(raw_forms, mng_ptr, block_str);
    std::string folderName = mng_ptr->prsm_para_ptr_->getOriDbName();
    folderName = folderName + "_index";
    //createIndexFiles<MassOnePtmFilterPtr>(block_str, folderName, filter_ptr);
}

void createIndexPtr(ProteoformPtrVec &raw_forms, int block_idx, DiagFilterMngPtr mng_ptr){
    std::string block_str = str_util::toString(block_idx);
    MassDiagFilterPtr filter_ptr = std::make_shared<MassDiagFilter>(raw_forms, mng_ptr, block_str);
    std::string folderName = mng_ptr->prsm_para_ptr_->getOriDbName();
    folderName = folderName + "_index";
    createIndexFiles<MassDiagFilterPtr>(block_str, folderName, filter_ptr);
}

void createIndexFiles(str::string block_idx, ZeroPtmFilterMngPtr mng_ptr){
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

std::function<void()> geneTask(int block_idx, DiagFilterMngPtr mng_ptr){
    return[block_idx, mng_ptr](){

    std::string db_block_file_name = mng_ptr->prsm_para_ptr_->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);    

    ProteoformPtrVec raw_forms 
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name, mng_ptr->prsm_para_ptr_->getFixModPtrVec());

    createIndexPtr(raw_forms, block_idx, mng_ptr);
    };
}
std::function<void()> geneTask(int block_idx, OnePtmFilterMngPtr mng_ptr){
    return[block_idx, mng_ptr](){

    std::string db_block_file_name = mng_ptr->prsm_para_ptr_->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);    

    ProteoformPtrVec raw_forms 
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name, mng_ptr->prsm_para_ptr_->getFixModPtrVec());

    createIndexPtr(raw_forms, block_idx, mng_ptr);
    };
}
*/
std::function<void()> geneTask(int block_idx, PrsmParaPtr prsm_para_ptr){
    return[block_idx, prsm_para_ptr](){

    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);    

    ProteoformPtrVec raw_forms 
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name, prsm_para_ptr->getFixModPtrVec());

    createIndexPtr(raw_forms, block_idx, prsm_para_ptr);
    };
}

void process(PrsmParaPtr prsm_para_ptr, int thread_num){

    std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
    DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

    int block_num = db_block_ptr_vec.size();

    SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);
    
    for (int i = 0; i < block_num; i++) {
        while (pool_ptr->getQueueSize() >= thread_num * 2) {
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
        pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), prsm_para_ptr));
  }
  pool_ptr->ShutDown();
}

}