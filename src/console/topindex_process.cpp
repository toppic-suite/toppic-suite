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
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "console/topindex_argument.hpp"
#include "console/topindex_process.hpp"

#include "common/base/base_data.hpp"
#include "common/util/file_util.hpp"

#include "seq/fasta_util.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

#include "filter/mng/zero_ptm_filter_mng.hpp"
#include "filter/mng/one_ptm_filter_mng.hpp"
#include "filter/mng/diag_filter_mng.hpp"
#include "filter/mng/topindex_file_name.hpp"
#include "filter/zeroindex/zero_ptm_index_processor.hpp"
#include "filter/oneindex/one_ptm_index_processor.hpp"
#include "filter/diagindex/diag_index_processor.hpp"

namespace toppic{

void TopIndexProcess(std::map<std::string, std::string> &arguments){
  try {
    Argument::outputArguments(std::cout, arguments);
    base_data::init();

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name = arguments["databaseFileName"];

    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    int thread_num = std::stoi(arguments["threadNumber"]);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);

    if (!file_util::exists(ori_db_file_name + "_idx")){//if _html folder was not created with topfd
      file_util::createFolder(ori_db_file_name + "_idx");
    }
    //prsm_para_ptr->setIndexDir(ori_db_file_name + "_idx"); //this gives seg fault

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    //create a folder for index files 
    // index file name
    TopIndexFileNamePtr file_name_ptr = std::make_shared<TopIndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);

    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);

    ZeroPtmFilterMngPtr zero_filter_mng_ptr
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, index_file_para, 
                                             thread_num, "");
    ZeroPtmIndexProcessorPtr zero_ptm_index_processor
        = std::make_shared<ZeroPtmIndexProcessor>(zero_filter_mng_ptr);

    zero_ptm_index_processor->process();
    zero_ptm_index_processor = nullptr;

    OnePtmFilterMngPtr one_ptm_filter_mng_ptr
        = std::make_shared<OnePtmFilterMng>(prsm_para_ptr, index_file_para, 
                                            "toppic_one_filter", thread_num);
    OnePtmIndexProcessorPtr one_ptm_index_processor
        = std::make_shared<OnePtmIndexProcessor>(one_ptm_filter_mng_ptr);
    one_ptm_index_processor->process();
    one_ptm_index_processor = nullptr;

    DiagFilterMngPtr diag_filter_mng_ptr
        = std::make_shared<DiagFilterMng>(prsm_para_ptr, index_file_para, 
                                          filter_result_num, thread_num, 
                                          "toppic_multi_filter");
    DiagIndexProcessorPtr diag_index_processor
        = std::make_shared<DiagIndexProcessor>(diag_filter_mng_ptr);
    diag_index_processor->process();
    diag_index_processor = nullptr;

    std::cout << "Deleting temporary files - started." << std::endl;

    std::string fa_base = file_util::absoluteName(ori_db_file_name);
    std::replace(fa_base.begin(), fa_base.end(), '\\', '/');
    file_util::cleanPrefix(ori_db_file_name, fa_base + "_");
    
    std::cout << "Deleting temporary files - finished." << std::endl; 

    std::cout << "TopIndex - finished." << std::endl;
    } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
}


    
}
