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

#include "util/file_util.hpp"
#include "util/str_util.hpp"
#include "util/logger.hpp"
#include "base/base_data.hpp"
#include "feature/deconv_process.hpp"
#include "feature/feature_detect.hpp"

namespace toppic {

int TopFDProcess(std::map<std::string, std::string> arguments) {
  try {
    time_t start = time(0);

    base_data::init(arguments["resourceDir"]);

    DeconvParaPtr para_ptr = std::make_shared<DeconvPara>(arguments);
    LOG_DEBUG("deconv para");
    DeconvProcess process(para_ptr);
    LOG_DEBUG("init process");
    process.process();
    FeatureDetect::process(para_ptr);

    time_t end = time(0);
    std::cout << "Runing time: "
        << str_util::toString(static_cast<int>(difftime(end, start)))
        << " seconds." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopFD finished." << std::endl << std::flush;
  return 0;
}

}  // namespace toppic

