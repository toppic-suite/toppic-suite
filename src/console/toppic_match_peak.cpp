//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include "base/version.hpp"
#include "base/base_data.hpp"
#include "prsm/prsm_coverage.hpp"

#include "console/toppic_argument.hpp"

namespace prot {

int outputMatchPeaks(std::map<std::string, std::string> arguments) {
  try {
    std::string exe_dir = arguments["executiveDir"];
    time_t start = time(0);
    char buf[50];
    arguments["start_time"] = std::string(ctime_r(&start, buf));
    Argument::outputArguments(std::cout, arguments);

    base_data::init(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    LOG_DEBUG("prsm para inited");

    std::cout << "Outputting coverage - started." << std::endl;
    PrsmCoveragePtr cov_ptr = std::make_shared<PrsmCoverage>(prsm_para_ptr, "CUTOFF_RESULT_SPEC", "COVERAGE");
    cov_ptr->processSingleCoverage();
    cov_ptr = nullptr;
    std::cout << "Outputting coverage - finished." << std::endl;


  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return EXIT_SUCCESS;
}

}

int main(int argc, char* argv[]) {
  // prot::log_level = 2;
  std::cout << std::setprecision(10);
  prot::Argument argu_processor;
  bool success = argu_processor.parse(argc, argv);
  if (!success) {
    return 1;
  }
  std::map<std::string, std::string> arguments = argu_processor.getArguments();
  prot::outputMatchPeaks(arguments);
}
