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
#include <map>
#include "console/topindex_argument.hpp"
#include "console/topindex_process.hpp"

#include "common/base/base_data.hpp"
#include "common/util/file_util.hpp"
#include "common/util/version.hpp"

#include "seq/fasta_util.hpp"
#include "seq/db_block.hpp"

#include "filter/mng/zero_ptm_filter_mng.hpp"
#include "filter/mng/one_ptm_filter_mng.hpp"
#include "filter/mng/diag_filter_mng.hpp"
#include "filter/mng/index_file_name.hpp"
#include "filter/index/zero_ptm_index.hpp"
#include "filter/index/one_ptm_index.hpp"
#include "filter/index/diag_index.hpp"

namespace toppic{
void TopIndexProcess(std::map<std::string, std::string> &arguments){
  try {
    std::cout << "TopIndex " << Version::getVersion() << std::endl;
    arguments["version"] = toppic::Version::getVersion();
    TopIndexArgument::outputArguments(std::cout, " ", arguments);
    base_data::init();

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name = arguments["databaseFileName"];

    int thread_num = std::stoi(arguments["threadNumber"]);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);

    if (!file_util::exists(ori_db_file_name + "_idx")){
      file_util::createFolder(ori_db_file_name + "_idx");
    }

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    // create a folder for index files 
    // index file name
    IndexFileNamePtr file_name_ptr = std::make_shared<IndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    int max_frag_len = std::stoi(arguments["maxFragmentLength"]);
    int min_block_num = std::stoi(arguments["minBlockNum"]);
    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, 
                             db_block_size, max_frag_len, min_block_num);

    ZeroPtmFilterMngPtr zero_filter_mng_ptr
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, index_file_para, 
                                             thread_num, "");
    zero_ptm_index::process(zero_filter_mng_ptr);

    OnePtmFilterMngPtr one_ptm_filter_mng_ptr
        = std::make_shared<OnePtmFilterMng>(prsm_para_ptr, index_file_para, 
                                            "toppic_one_filter", thread_num);
    one_ptm_index::process(one_ptm_filter_mng_ptr);

    DiagFilterMngPtr diag_filter_mng_ptr
        = std::make_shared<DiagFilterMng>(prsm_para_ptr, index_file_para, 
                                          filter_result_num, thread_num, 
                                          "toppic_multi_filter");
    diag_index::process(diag_filter_mng_ptr);

    base_data::release();

    std::cout << "TopIndex - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
}

}
