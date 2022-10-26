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

#include <algorithm>
#include <iostream>
#include <sstream>

#include "common/util/logger.hpp"
#include "gui/util/run_exe.hpp"


namespace toppic {

namespace run_exe {

/*function for topfd*/ 
std::string geneTopfdCommand(TopfdParaPtr para_ptr, 
                             const std::vector<std::string> spec_file_lst, 
                             std::string app_name) {

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = para_ptr->getExeDir() + "\\" + app_name + ".exe ";
#else
  std::string exe_path = para_ptr->getExeDir() + "/" + app_name + " ";
#endif

  std::string command = exe_path;
  std::stringstream oss;
  oss << "-a " << para_ptr->getActivation() << " ";
  oss << "-c " << para_ptr->getMaxCharge() << " ";
  oss << "-m " << para_ptr->getMaxMass() << " ";
  oss << "-t " << para_ptr->getMzError() << " ";
  oss << "-r " << para_ptr->getMsOneSnRatio() << " ";
  oss << "-s " << para_ptr->getMsTwoSnRatio() << " ";
  oss <<  "-w " << para_ptr->getPrecWindow() << " ";
  command = command + oss.str();
  if (para_ptr->isUseEnvCnn()) {
    command = command + "-n ";
  }
  if (para_ptr->isMissingLevelOne()) {
    command = command + "-o ";
  }
  command = command + "-u " + std::to_string(para_ptr->getThreadNum()) + " ";
  if (!para_ptr->isGeneHtmlFolder()) {
    command = command + "-g ";
  }
  if (!para_ptr->isDoFinalFiltering()) {
    command = command + "-d ";
  }
  for (size_t i = 0; i < spec_file_lst.size(); i++) {
    command = command + spec_file_lst[i] + " ";
  }
  return command;
}

std::map<std::string, std::string> topindex_para {
  {"fixedMod", "-f"},
    {"allowProtMod", "-n"},
    {"searchType", "-d"},
    {"threadNumber", "-u"},
    {"massErrorTolerance", "-e"}
};


/*function for topindex*/
std::string geneTopIndexCommand(std::map<std::string, std::string> arguments_, 
                                std::string app_name) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + app_name + " ";
  #endif

  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    //if one of the topindex parameters
    if (topindex_para.find(it->first) != topindex_para.end()) { 
      //some parameters require extra processing
      if (it->first == "fixedMod" && it->second == "") {
        continue; //don't add -f
      }
      else if (it->first == "searchType") {
        if (it->second != "TARGET") {
          command = command + "-d ";
        }
        continue;
      }
      command = command + topindex_para[it->first] + " " + it->second + " ";
    }
  }
  command = command + arguments_["oriDatabaseFileName"] + " ";
  return command;
};

std::map<std::string, std::string> toppic_para {
  {"activation", "-a "},
    {"fixedMod", "-f "},
    {"allowProtMod", "-n "},
    {"searchType", "-d "},
    {"threadNumber", "-u "},
    {"proteoformErrorTolerance", "-p "},
    {"maxPtmMass", "-M "},
    {"minPtmMass", "-m "},
    {"cutoffSpectralType", "-t "},
    {"cutoffSpectralValue", "-v "},
    {"cutoffProteoformType", "-T "},
    {"cutoffProteoformValue", "-V "},
    {"ptmNumber", "-s "},
    {"useFeatureFile", "-x "},
    {"keepTempFiles", "-k "},
    {"keepDecoyResults", "-K "},
    {"geneHTMLFolder", "-g "},
    {"combinedOutputName", "-c "},
    {"massErrorTolerance", "-e "},
    {"useLookupTable", "-l "},
    {"groupSpectrumNumber", "-r "},
    {"residueModFileName", "-i "},
    {"localThreshold", "-H "}
};

std::string geneToppicCommand(std::map<std::string, std::string> arguments_, 
                              std::vector<std::string> spec_file_lst_, 
                              std::string app_name) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + app_name + " ";
  #endif

  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    //if one of the toppic parameters
    if (toppic_para.find(it->first) != toppic_para.end()) { 
      //skip some paramters based on parameter values
      if (it->first == "fixedMod" && it->second == "") {
        continue;
      }
      else if (it->first == "combinedOutputName" && it->second == "") {
        continue;
      }
      else if (it->first == "useFeatureFile") {
        if (it->second == "false") {
          command = command + toppic_para[it->first] + " ";
        }
      }
      else if (it->first == "searchType") {
        if (it->second != "TARGET") {
          command = command + toppic_para[it->first] + " ";
        }
      }
      else if (it->first == "keepTempFiles" || it->first == "keepDecoyResults") {
        if (it->second == "true") {
          command = command + toppic_para[it->first] + " ";
        }
      }
      else if (it->first == "geneHTMLFolder" ) {//for geneHTML folder, the argument should be added when the value is false
        if (it->second != "true") {
          command = command + toppic_para[it->first] + " ";
        }
      }
      else if (it->first == "residueModFileName" && it->second == "") {
        continue; //don't add -i
      }
      else if (it->first == "useLookupTable") {
        if (it->second == "true") {
          command = command + toppic_para[it->first] + " ";
        }
      }
      else{
        command = command + toppic_para[it->first] + it->second + " ";
      }
    }
    else {//parameter is not found anywhere
      LOG_DEBUG("Parameter " << it->first << " from " << app_name << " was not found in any apps!");
    }
  }  
  command = command + arguments_["oriDatabaseFileName"] + " ";

  for (size_t i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
};

std::map<std::string, std::string> topmg_para {
  {"activation", "-a "},
    {"fixedMod", "-f "},
    {"allowProtMod", "-n "},
    {"searchType", "-d "},
    {"threadNumber", "-u "},
    {"proteoformErrorTolerance", "-p "},
    {"maxPtmMass", "-M "},
    {"cutoffSpectralType", "-t "},
    {"cutoffSpectralValue", "-v "},
    {"cutoffProteoformType", "-T "},
    {"cutoffProteoformValue", "-V "},
    {"ptmNumber", "-s "},
    {"useFeatureFile", "-x "},
    {"keepTempFiles", "-k "},
    {"keepDecoyResults", "-K "},
    {"geneHTMLFolder", "-g "},
    {"combinedOutputName", "-c "},
    {"massErrorTolerance", "-e "},
    {"useAsfDiag", "-D "},
    {"varModFileName", "-i "},
    {"varPtmNumber", "-P "},
    {"wholeProteinOnly", "-w "},
    {"proteoGraphGap", "-j "},
    {"varPtmNumInGap", "-G "}
};

std::string geneTopmgCommand(std::map<std::string, std::string> arguments_, 
                             std::vector<std::string> spec_file_lst_, 
                             std::string app_name) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + app_name + " ";
  #endif

  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    if (topmg_para.find(it->first) != topmg_para.end()) { //if one of the common parameters
      //skip some paramters based on parameter values
      if (it->first == "fixedMod" && it->second == "") {
        continue;
      }
      else if (it->first == "combinedOutputName" && it->second == "") {
        continue;
      }
      else if (it->first == "useFeatureFile") {
        if (it->second == "false") {
          command = command + topmg_para[it->first] + " ";
        }
      }
      else if (it->first == "searchType") {
        if (it->second != "TARGET") {
          command = command + topmg_para[it->first] + " ";
        }
      }
      else if (it->first == "keepTempFiles" || it->first == "keepDecoyResults") {
        if (it->second == "true") {
          command = command + topmg_para[it->first] + " ";
        }
      }
      else if (it->first == "geneHTMLFolder" ) {//for geneHTML folder, the argument should be added when the value is false
        if (it->second != "true") {
          command = command + topmg_para[it->first] + " ";
        }
      }
      //some parameters require extra processing
      else if (it->first == "useAsfDiag" || it->first == "wholeProteinOnly") {
        if (it->second == "true") {
          command = command + topmg_para[it->first] + " ";
        }
      }
      else {
        command = command + topmg_para[it->first] + it->second + " ";
      }
    }
    else {//parameter is not found anywhere
      LOG_DEBUG("Parameter " << it->first << " from " << app_name << " was not found in any apps!");
    }
  }  
  command = command + arguments_["oriDatabaseFileName"] + " ";

  for (size_t i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
};

std::map<std::string, std::string> topdiff_para {
  {"errorTolerance", "-e "},
    {"mergedOutputFileName", "-o "},
    {"toolName", "-t "},
};


// function for topdiff
std::string geneTopDiffCommand(std::map<std::string, std::string> arguments_, 
                               std::vector<std::string> spec_file_lst_, 
                               std::string app_name) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + app_name + ".exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + app_name + " ";
  #endif

  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    if (app_name == "topdiff" && topdiff_para.find(it->first) != topdiff_para.end()) {
      command = command + topdiff_para[it->first] + it->second + " ";
    }
    else {//parameter is not found anywhere
      LOG_DEBUG("Parameter " << it->first << " from " << app_name << " was not found in any apps!");
    }
  }  
  for (size_t i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
};

void startJob() {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  HANDLE hJob = CreateJobObject( NULL, NULL );
  if ( hJob == NULL ) {
    std::cout << "CreateJobObject failed: error " << GetLastError() << std::endl;
    return;
  }

  JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = { 0 };
  jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
  bool bSuccess = SetInformationJobObject( hJob, JobObjectExtendedLimitInformation, &jeli, sizeof( jeli ) );
  if ( bSuccess == 0 ) {
    std::cout << "SetInformationJobObject failed: error " << GetLastError() << std::endl;
    return;
  }
  bSuccess = AssignProcessToJobObject(hJob, GetCurrentProcess());
  if ( bSuccess == 0) {
     std::cout << "AssignProcessToJobObject failed!" << std::endl;
     return;
  }
  #endif
}

void run(std::string command) {
  std::cout << command << std::endl;
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
  bSuccess = CreateProcess(NULL, command.c_str(), NULL, NULL, 
                           TRUE, CREATE_NO_WINDOW, NULL, NULL, 
                           &siStartInfo, &piProcInfo);  // receives PROCESS_INFORMATION 

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
    char buf[4096]; 
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
      LOG_ERROR("Error occured when opening pipe");
    }
    while (fgets(buf, 4096, pipe) != NULL) {
        std::cout << buf;
    }
    pclose(pipe);
  #endif
};

}

}
