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

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/base/base_data.hpp"
#include "seq/fasta_util.hpp"
#include "merge/feature_sample_merge.hpp"
#include "console/topdiff_argument.hpp"
#include "console/topdiff_process.hpp"

namespace toppic {

int topDiffProcess(std::map<std::string, std::string> &arguments,
                    std::vector<std::string> &input_file_list) {
  try{
    Argument::outputArguments(std::cout, arguments);
    base_data::init();

    std::string ori_db_file_name = arguments["databaseFileName"];
    std::string db_file_name = ori_db_file_name + "_target";
    fasta_util::dbSimplePreprocess(ori_db_file_name, db_file_name);

    std::string base_path = file_util::absoluteDir(input_file_list[0]);
    std::string output_file_name = base_path + file_util::getFileSeparator() 
        + arguments["mergedOutputFileName"];
    LOG_DEBUG("Output file name " << output_file_name);
    std::string fixed_mod = arguments["fixedMod"];

    double error_tole = std::stod(arguments["errorTolerance"]);
    std::string tool_name = arguments["toolName"];

    std::cout << "TopDiff - started." << std::endl;
    FeatureSampleMergePtr merge_ptr 
        = std::make_shared<FeatureSampleMerge>(input_file_list,
                                                output_file_name,
                                                db_file_name, 
                                                fixed_mod, 
                                                tool_name,
                                                error_tole);
    merge_ptr->process();
    merge_ptr = nullptr;
    
    std::cout << "Deleting temporary files - started." << std::endl;

    std::string fa_base = file_util::absoluteName(ori_db_file_name);
    std::replace(fa_base.begin(), fa_base.end(), '\\', '/');
    file_util::cleanPrefix(ori_db_file_name, fa_base + "_");
        
    std::cout << "Deleting temporary files - finished." << std::endl; 

    std::cout << "TopDiff - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}

}
