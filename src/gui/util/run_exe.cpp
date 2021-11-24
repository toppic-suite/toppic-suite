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

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#include <windows.h> 
#include <tchar.h>
#include <strsafe.h>
#else
#include <array>
#endif

#include "common/util/logger.hpp"
#include "gui/util/run_exe.h"
#include <algorithm>
#include <iostream>
#include <stdio.h> 

namespace toppic {
/*function for topindex*/
std::string RunExe::geneCommand(std::map<std::string, std::string> arguments_, std::string app_name) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + app_name + " ";
  #endif

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
  command = command + arguments_["oriDatabaseFileName"] + " ";
  return command;
};
/*function for topfd*/ 
std::string RunExe::geneCommand(TopfdParaPtr para_ptr, std::vector<std::string> spec_file_lst_, std::string app_name) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = para_ptr->exe_dir_ + "\\" + app_name + ".exe ";
  #else
  std::string exe_path = para_ptr->exe_dir_ + "/" + app_name + " ";
  #endif

  std::string command = exe_path;
  command = command + "-c " + std::to_string(para_ptr->max_charge_) + " ";
  command = command + "-m " + std::to_string(para_ptr->max_mass_) + " ";
  command = command + "-t " + std::to_string(para_ptr->mz_error_) + " ";
  command = command + "-r " + std::to_string(para_ptr->ms_one_sn_ratio_) + " ";
  command = command + "-s " + std::to_string(para_ptr->ms_two_sn_ratio_) + " ";
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
  for (size_t i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
}

/*function for toppic, topmg, topmerge, topdiff*/
std::string RunExe::geneCommand(std::map<std::string, std::string> arguments_, std::vector<std::string> spec_file_lst_, std::string app_name) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + app_name + " ";
  #endif

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
      else if (it->first == "keepTempFiles" || it->first == "keepDecoyResults") {
        if (it->second == "true") {
          command = command + common_para[it->first] + " ";
        }
      }
      else if (it->first == "geneHTMLFolder" ) {//for geneHTML folder, the argument should be added when the value is false
        if (it->second != "true") {
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
    else if (app_name == "topmerge") {
      //some parameters require extra processing
      if (it->first == "residueModFileName" && it->second == "") continue; //don't add -i
      else if (it->first == "useLookupTable") {
        if (it->second == "true") {
          command = command + topmerge_para[it->first] + " ";
        }
      }
      else {
        command = command + topmerge_para[it->first] + it->second + " ";
      }
    }
    else {//parameter is not found anywhere
      LOG_ERROR("Parameter " << it->first << " from " << app_name << " was not found in any apps!");
      return "";
    }
  }  
  command = command + arguments_["oriDatabaseFileName"] + " ";

  for (size_t i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
};
void RunExe::run(std::string command) {
  //std::cout << command << std::endl;
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  HANDLE g_hChildStd_IN_Rd = NULL;
  HANDLE g_hChildStd_IN_Wr = NULL;
  HANDLE g_hChildStd_OUT_Rd = NULL;
  HANDLE g_hChildStd_OUT_Wr = NULL;

  SECURITY_ATTRIBUTES saAttr; 

  //Set the bInheritHandle flag so pipe handles are inherited. 
  saAttr.nLength = sizeof(SECURITY_ATTRIBUTES); 
  saAttr.bInheritHandle = TRUE; 
  saAttr.lpSecurityDescriptor = NULL; 

  CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &saAttr, 0);
  CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &saAttr, 0);

  SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0);
  SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0);

  PROCESS_INFORMATION piProcInfo; 
  STARTUPINFO siStartInfo;
  BOOL bSuccess = FALSE; 
 
  ZeroMemory( &piProcInfo, sizeof(PROCESS_INFORMATION) );
  ZeroMemory( &siStartInfo, sizeof(STARTUPINFO));

  siStartInfo.cb = sizeof(STARTUPINFO); 
  siStartInfo.hStdError = g_hChildStd_OUT_Wr;
  siStartInfo.hStdOutput = g_hChildStd_OUT_Wr;
  siStartInfo.hStdInput = g_hChildStd_IN_Rd;
  siStartInfo.dwFlags |= STARTF_USESTDHANDLES;
 
  //Create the child process. 
  bSuccess = CreateProcess(NULL, command.c_str(), NULL, NULL, TRUE, CREATE_NO_WINDOW, NULL, NULL, &siStartInfo, &piProcInfo);  // receives PROCESS_INFORMATION 

  if (!bSuccess ) {
    std::cout << "error occured" << std::endl;
    std::cout << "command: " << command << std::endl;
    return;
  }
  else {
    CloseHandle(piProcInfo.hProcess);
    CloseHandle(piProcInfo.hThread);
    CloseHandle(g_hChildStd_OUT_Wr);
    CloseHandle(g_hChildStd_IN_Rd);
  }   
  DWORD dwRead, dwWritten; 
  CHAR buf[4096]; 
  BOOL readSuccess = FALSE;
  HANDLE hParentStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
  bSuccess = ReadFile( g_hChildStd_OUT_Rd, buf, 4096, &dwRead, NULL);
  std::ios::sync_with_stdio(true);
  while (bSuccess == TRUE) {
    buf[dwRead] = '\0';
    OutputDebugStringA(buf);
    std::cout << buf;
    bSuccess = ReadFile(g_hChildStd_OUT_Rd, buf, 1024, &dwRead, NULL);
  }
  #else
    //std::array<char, 128> buffer;
    char buf[4096]; 
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
      std::cout << "error occured when opening pipe" << std::endl;
    }
    while (fgets(buf, 4096, pipe) != NULL) {
        std::cout << buf;
    }
    pclose(pipe);
  #endif
};
}
