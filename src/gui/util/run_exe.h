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

#ifndef GUI_RUN_EXE_H
#define GUI_RUN_EXE_H

#include "topfd/common/topfd_para.hpp"
#include <vector>
#include <map>
#include <string>

namespace toppic {
class RunExe {
 public:
  std::map<std::string, std::string> common_para {
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
    {"combinedOutputName", "-c "}
    //{"oriDatabaseFileName", ""}
  };
  std::vector<std::string> skip_para {//parameters to skip
    "executiveDir", "resourceDir", "databaseBlockSize","filteringResultNumber",
    "groupSpectrumNumber","maxFragmentLength","numOfTopPrsms","skipList",
    "databaseFileName","oriDatabaseFileName","useGf"
  };
  std::map<std::string, std::string> topindex_para {
    {"massErrorTolerance", "-e "}
  };

  std::map<std::string, std::string> toppic_para {
    {"massErrorTolerance", "-e "},
    {"useLookupTable", "-l "},
    {"groupSpectrumNumber", "-r "},
    {"residueModFileName", "-i "},
    {"localThreshold", "-H "}
  };

  std::map<std::string, std::string> topmg_para {
    {"massErrorTolerance", "-e "},
    {"useAsfDiag", "-D "},
    {"varModFileName", "-i "},
    {"varPtmNumber", "-P "},
    {"wholeProteinOnly", "-w "},
    {"proteoGraphGap", "-j "},
    {"varPtmNumInGap", "-G "}
  };

  std::map<std::string, std::string> topdiff_para {
    {"errorTolerance", "-e "},
    {"mergedOutputFileName", "-o "},
    {"toolName", "-t "},
  };

  std::map<std::string, std::string> topmerge_para {
    {"residueModFileName", "-i "},
    {"localThreshold", "-H "},
    {"massErrorTolerance", "-e "}
  };

  std::string geneCommand(std::map<std::string, std::string> arguments_, std::string app_name);
  std::string geneCommand(TopfdParaPtr para_ptr, std::vector<std::string> spec_file_lst_, std::string app_name);
  std::string geneCommand(std::map<std::string, std::string> arguments_, std::vector<std::string> spec_file_lst_, std::string app_name);
  void run(std::string command); 
};
}
#endif