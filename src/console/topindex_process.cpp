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
#include <ctime>

#include "console/topindex_argument.hpp"
#include "console/topindex_process.hpp"

#include "common/base/base_data.hpp"
#include "common/util/file_util.hpp"

#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "filter/diag/diag_filter_mng.hpp"
#include "filter/diag/diag_filter_processor.hpp"

#include "seq/fasta_util.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic{

void TopIndexProcess(std::map<std::string, std::string> &arguments){
    Argument::outputArguments(std::cout, arguments);

    base_data::init();

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name = arguments["databaseFileName"];

    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    int thread_num = std::stoi(arguments["threadNumber"]);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);

    file_util::createFolder(ori_db_file_name + "_idx");
    //prsm_para_ptr->setIndexDir(ori_db_file_name + "_idx"); //this gives seg fault

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    
    //create a folder for index files 

    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);

    ZeroPtmFilterMngPtr zero_filter_mng_ptr
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, thread_num, "");
    ZeroPtmFilterProcessorPtr zero_filter_processor
        = std::make_shared<ZeroPtmFilterProcessor>(zero_filter_mng_ptr);

    zero_filter_processor->index_process();
    zero_filter_processor = nullptr;

    OnePtmFilterMngPtr one_ptm_filter_mng_ptr
        = std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "toppic_one_filter", thread_num);
    OnePtmFilterProcessorPtr one_filter_processor
        = std::make_shared<OnePtmFilterProcessor>(one_ptm_filter_mng_ptr);
    one_filter_processor->index_process();
    one_filter_processor = nullptr;

    DiagFilterMngPtr diag_filter_mng_ptr
        = std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num,
                                            thread_num, "toppic_multi_filter");
    DiagFilterProcessorPtr diag_filter_processor
        = std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
    diag_filter_processor->index_process();
    diag_filter_processor = nullptr;
    
}

}