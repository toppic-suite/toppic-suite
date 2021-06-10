//Copyright (c) 2014 - 2021, The Trustees of Indiana University.
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
#include "gui/util/run_program.h"
#include <vector>

namespace toppic {

std::string RunProgram::geneCommand(std::map<std::string, std::string> arguments_, std::vector<std::string> spec_file_lst_, std::string app_name) {
  //have each argument converted to command line parameter
  //create a string for a command at the end
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    if (it->first == "executiveDir" || it->first == "resourceDir") continue;
    
    else if (common_para.find(it->first) != common_para.end()) { //if one of the common parameters
      //some parameters require extra processing
      if (it->first == "fixedMod" && it->second == "") continue; //don't add -f
      if (it->first == "searchType" && it->second == "TARGET") continue; // don't add -d
      command = command + common_para[it->first] + it->second + " ";
    }
    else if (app_name == "topfd") {}// this needs separate processing
    else if (app_name == "toppic" && toppic_para.find(it->first) != toppic_para.end()) {
      command = command + toppic_para[it->first] + it->second + " ";
    }
    else if (app_name == "topindex" && topindex_para.find(it->first) != topindex_para.end()) {
      command = command + topindex_para[it->first] + it->second + " ";
    }
    else if (app_name == "topmg" && topmg_para.find(it->first) != topmg_para.end()) {
      command = command + topmg_para[it->first] + it->second + " ";
    }
    else if (app_name == "topdiff" && topdiff_para.find(it->first) != topdiff_para.end()) {
      command = command + topdiff_para[it->first] + it->second + " ";
    }
    else {//parameter is not found anywhere
      LOG_ERROR("Parameter " << it->first << " from " << app_name << " was not found in any apps!");
      return "";
    }
  }
  if (app_name == "topfd") {}
  else if (app_name == "toppic") {}
  else if (app_name == "topindex") {}
  else if (app_name == "topmg") {}
  else if (app_name == "topdiff") {
    for (int i = 0; i < spec_file_lst_.size(); i++) {
      command = command + spec_file_lst_[i] + " ";
    }
  }

  return command;
};
void RunProgram::run(std::string command) {
  if (command == "") {
    return;
  }
  std::cout << "cmd: " << command << std::endl;
  char buffer[500];
  FILE* pipe = _popen(command.c_str(), "r");
  if (!pipe) {
      LOG_ERROR("popen failed!");
  }
  while (!feof(pipe)) {
    if (fgets(buffer, 128, pipe) != NULL) {
        std::cout << buffer;
    }
  }
  _pclose(pipe);
};
}
