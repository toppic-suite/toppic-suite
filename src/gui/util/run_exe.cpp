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
#include "gui/util/run_exe.h"
#include <algorithm>
#include <iostream>

namespace toppic {
/*function for topindex*/
std::string RunExe::geneCommand(std::map<std::string, std::string> arguments_, std::string app_name) {
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    if (std::count(skip_para.begin(), skip_para.end(), it->first)) continue;
    else if (common_para.find(it->first) != common_para.end()) { //if one of the common parameters
      //some parameters require extra processing
      if (it->first == "fixedMod" && it->second == "") continue; //don't add -f
      else if (it->first == "activation") continue; //for topindex, activation para should be skipped
      else if (it->first == "searchType") {
        if (it->second != "TARGET") {
          command = command + "-d ";
        }
        continue;
      }
      command = command + common_para[it->first] + it->second + " ";
    }
    else if (app_name == "topindex" && topindex_para.find(it->first) != topindex_para.end()) {
      command = command + topindex_para[it->first] + it->second + " ";
    }
    else {//parameter is not found anywhere
      LOG_ERROR("Parameter " << it->first << " from " << app_name << " was not found in any apps!");
      return "";
    }
  }
  return command;
};
/*function for topfd*/ 
std::string RunExe::geneCommand(TopfdParaPtr para_ptr, std::vector<std::string> spec_file_lst_, std::string app_name) {
  std::string exe_path = para_ptr->exe_dir_ + "\\" + app_name + ".exe ";
  std::string command = exe_path;

  command = command + "-c " + std::to_string(para_ptr->max_charge_) + " ";
  command = command + "-m " + std::to_string(para_ptr->max_mass_) + " ";
  command = command + "-e " + std::to_string(para_ptr->mz_error_) + " ";
  command = command + "-r " + std::to_string(para_ptr->ms_one_sn_ratio_) + " ";
  command = command + "-t " + std::to_string(para_ptr->ms_two_sn_ratio_) + " ";
  command = command + "-w " + std::to_string(para_ptr->prec_window_) + " ";
  command = command + "-u " + std::to_string(para_ptr->thread_number_) + " ";
  command = command + "-a " + para_ptr->activation_ + " ";
  
  if (para_ptr->missing_level_one_) {
    command = command + "-o ";
  }
  if (!para_ptr->gene_html_folder_) {
    command = command + "-g ";
  }
  if (para_ptr->use_env_cnn_) {
    command = command + "-n ";
  }
  if (!para_ptr->do_final_filtering_) {
    command = command + "-d ";
  }
  for (int i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
}

/*function for toppic, topmg, topdiff*/
std::string RunExe::geneCommand(std::map<std::string, std::string> arguments_, std::vector<std::string> spec_file_lst_, std::string app_name) {
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    if (std::count(skip_para.begin(), skip_para.end(), it->first)) continue;
    else if (common_para.find(it->first) != common_para.end()) { //if one of the common parameters
      //skip some paramters based on parameter values
      if (it->first == "fixedMod" && it->second == "") continue;
      else if (it->first == "combinedOutputName" && it->second == "") continue;
      else if (it->first == "useFeatureFile") {
        if (it->second == "false") {
          command = command + common_para[it->first] + " ";
        }
      }
      else if (it->first == "searchType") {
        if (it->second != "TARGET") {
          command = command + common_para[it->first] + " ";
        }
      }
      else if (it->first == "keepTempFiles" || it->first == "geneHTMLFolder") {
        if (it->second == "true") {
          command = command + common_para[it->first] + " ";
        }
      }
      else{
        command = command + common_para[it->first] + it->second + " ";
      }
    }
    else if (app_name == "toppic" && toppic_para.find(it->first) != toppic_para.end()) {
      //some parameters require extra processing
      if (it->first == "residueModFileName" && it->second == "") continue; //don't add -i
      else if (it->first == "useLookupTable") {
        if (it->second == "true") {
          command = command + toppic_para[it->first] + " ";
        }
      }
      else {
        command = command + toppic_para[it->first] + it->second + " ";
      }
    }
    else if (app_name == "topmg" && topmg_para.find(it->first) != topmg_para.end()) {
      //some parameters require extra processing
      if (it->first == "useAsfDiag" || it->first == "wholeProteinOnly") {
        if (it->second == "true") {
          command = command + topmg_para[it->first] + " ";
        }
      }
      else {
        command = command + topmg_para[it->first] + it->second + " ";
      }
    }
    else if (app_name == "topdiff" && topdiff_para.find(it->first) != topdiff_para.end()) {
      command = command + topdiff_para[it->first] + it->second + " ";
    }
    else {//parameter is not found anywhere
      LOG_ERROR("Parameter " << it->first << " from " << app_name << " was not found in any apps!");
      return "";
    }
  }  
  command = command + arguments_["oriDatabaseFileName"] + " ";

  for (int i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
};
void RunExe::run(std::string command) {
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
