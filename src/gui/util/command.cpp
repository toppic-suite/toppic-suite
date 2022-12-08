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

#include <iostream>
#include <sstream>

#include "common/util/logger.hpp"
#include "gui/util/command.hpp"


namespace toppic {

namespace command {

/*function for topfd*/ 
std::string geneTopfdCommand(TopfdParaPtr para_ptr, 
                             const std::vector<std::string> spec_file_lst) { 

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = para_ptr->getExeDir() + "\\" + "topfd.exe ";
#else
  std::string exe_path = para_ptr->getExeDir() + "/" + "topfd ";
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
std::string geneTopIndexCommand(std::map<std::string, 
                                std::string> arguments_) {
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + "topindex.exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + "topindex ";
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
    {"maxShiftMass", "-M "},
    {"minShiftMass", "-m "},
    {"variablePtmNum", "-S "},
    {"variablePtmFileName", "-b "},
    {"useApproxSpectra", "-A "},
    {"cutoffSpectralType", "-t "},
    {"cutoffSpectralValue", "-v "},
    {"cutoffProteoformType", "-T "},
    {"cutoffProteoformValue", "-V "},
    {"shiftNumber", "-s "},
    {"useFeatureFile", "-x "},
    {"keepTempFiles", "-k "},
    {"keepDecoyResults", "-K "},
    {"geneHTMLFolder", "-g "},
    {"combinedOutputName", "-c "},
    {"massErrorTolerance", "-e "},
    {"useLookupTable", "-l "},
    {"groupSpectrumNumber", "-r "},
    {"localPtmFileName", "-B "},
    {"localThreshold", "-H "}
};

std::string geneToppicCommand(std::map<std::string, std::string> arguments_, 
                              std::vector<std::string> spec_file_lst_) { 
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + "toppic.exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + "toppic ";
  #endif

  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    //if one of the toppic parameters
    if (toppic_para.find(it->first) != toppic_para.end()) { 
      //skip some paramters based on parameter values
      LOG_DEBUG(it->first << " " << it->second);
      if (it->first == "fixedMod" && it->second == "") {
        continue;
      }
      else if (it->first == "combinedOutputName" && it->second == "") {
        continue;
      }
      else if (it->first == "useFeatureFile") {
        if (it->second == "false") {
          command = command + toppic_para[it->first];
        }
      }
      else if (it->first == "searchType") {
        if (it->second != "TARGET") {
          command = command + toppic_para[it->first];
        }
      }
      else if (it->first == "keepTempFiles" || it->first == "keepDecoyResults") {
        if (it->second == "true") {
          command = command + toppic_para[it->first];
        }
      }
      else if (it->first == "geneHTMLFolder" ) {//for geneHTML folder, the argument should be added when the value is false
        if (it->second != "true") {
          command = command + toppic_para[it->first];
        }
      }
      else if (it->first == "useApproxSpectra") {
        if (it->second == "true") {
          command = command + toppic_para[it->first];
        }
      }
      else if (it->first == "variablePtmFileName" && it->second == "") {
        continue; //don't add -b
      }
      else if (it->first == "localPtmFileName" && it->second == "") {
        continue; //don't add -B
      }
      else if (it->first == "useLookupTable") {
        if (it->second == "true") {
          command = command + toppic_para[it->first];
        }
      }
      else{
        command = command + toppic_para[it->first] + it->second + " ";
      }
    }
    else {//parameter is not found anywhere
      LOG_DEBUG("Parameter " << it->first << " from toppic was not found in any apps!");
    }
  }  
  command = command + arguments_["oriDatabaseFileName"] + " ";

  for (size_t i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  LOG_DEBUG(command);
  return command;
};

std::map<std::string, std::string> topmg_para {
  {"activation", "-a "},
    {"fixedMod", "-f "},
    {"allowProtMod", "-n "},
    {"searchType", "-d "},
    {"threadNumber", "-u "},
    {"proteoformErrorTolerance", "-p "},
    {"maxShiftMass", "-M "},
    {"cutoffSpectralType", "-t "},
    {"cutoffSpectralValue", "-v "},
    {"cutoffProteoformType", "-T "},
    {"cutoffProteoformValue", "-V "},
    {"shiftNumber", "-s "},
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
                             std::vector<std::string> spec_file_lst_) { 
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + "topmg.exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + "topmg ";
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
          command = command + topmg_para[it->first];
        }
      }
      else if (it->first == "searchType") {
        if (it->second != "TARGET") {
          command = command + topmg_para[it->first];
        }
      }
      else if (it->first == "keepTempFiles" || it->first == "keepDecoyResults") {
        if (it->second == "true") {
          command = command + topmg_para[it->first];
        }
      }
      else if (it->first == "geneHTMLFolder" ) {//for geneHTML folder, the argument should be added when the value is false
        if (it->second != "true") {
          command = command + topmg_para[it->first];
        }
      }
      //some parameters require extra processing
      else if (it->first == "useAsfDiag" || it->first == "wholeProteinOnly") {
        if (it->second == "true") {
          command = command + topmg_para[it->first];
        }
      }
      else {
        command = command + topmg_para[it->first] + it->second + " ";
      }
    }
    else {//parameter is not found anywhere
      LOG_DEBUG("Parameter " << it->first << " from topmg was not found in any apps!");
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
                               std::vector<std::string> spec_file_lst_) { 
  #if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::string exe_path = arguments_["executiveDir"] + "\\" + "topdiff.exe ";
  #else
  std::string exe_path = arguments_["executiveDir"] + "/" + "topdiff ";
  #endif

  std::string command = exe_path;

  for (std::map<std::string, std::string>::iterator it = arguments_.begin(); it != arguments_.end(); ++it) {
    if (topdiff_para.find(it->first) != topdiff_para.end()) {
      command = command + topdiff_para[it->first] + it->second + " ";
    }
    else {//parameter is not found anywhere
      LOG_DEBUG("Parameter " << it->first << " from topdiff was not found in any apps!");
    }
  }  
  for (size_t i = 0; i < spec_file_lst_.size(); i++) {
    command = command + spec_file_lst_[i] + " ";
  }
  return command;
};

}

}
