//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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
#include "common/util/console_util.hpp"

namespace toppic {

namespace console_util {

void copyTopMSV(std::map<std::string, std::string> &arguments) {
  std::string spectrum_file_name = arguments["spectrumFileName"];
  std::string base_name = file_util::basename(spectrum_file_name);
  std::string base_name_short = base_name.substr(0, base_name.length() - 4);
  std::string topmsv_dir = base_name_short + "_html" +  file_util::getFileSeparator() + "topmsv";
  if (file_util::exists(topmsv_dir)) {
    LOG_WARN("The TopMSV directory " << topmsv_dir << " exists!");
    file_util::delDir(topmsv_dir);
  }
  if (!file_util::exists(base_name_short + "_html")){//if _html folder was not created with topfd
    file_util::createFolder(base_name_short + "_html");
  }
  std::string resource_dir = arguments["resourceDir"];
  // copy resources 
  std::string from_path(resource_dir + file_util::getFileSeparator() + "topmsv");
  file_util::copyDir(from_path, topmsv_dir);
}

}  // namespace console_util

}  // namespace toppic

