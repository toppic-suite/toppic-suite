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

#include <map>
#include <iostream>
#include <string>
#include <algorithm>

#include "console/topmerge_argument.hpp"
#include "console/topmerge_process.hpp"

using namespace toppic;

int main(int argc, char* argv[]) {
  //toppic::log_level = 2;
  LOG_DEBUG("Parsing start!");
  Argument argu_processor;
  bool success = argu_processor.parse(argc, argv);
  
  LOG_DEBUG("Parsing success!");

  if (!success) {
    return 1;
  }

  std::map<std::string, std::string> arguments = argu_processor.getArguments();

  std::vector<std::string> proteo_file_list = argu_processor.getProteoformFileList();

  topMergeProcess(arguments, proteo_file_list);

  return 0;
}
