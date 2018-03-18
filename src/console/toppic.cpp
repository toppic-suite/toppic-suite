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
#include <ctime>

#include "base/version.hpp"
#include "base/file_util.hpp"
#include "base/string_util.hpp"

#include "spec/msalign_util.hpp"
#include "spec/feature_util.hpp"

#include "prsm/prsm_util.hpp"

#include "console/toppic_argument.hpp"
#include "console/toppic_process.hpp"

int main(int argc, char* argv[]) {
  //prot::log_level = 2;
  std::cout << std::setprecision(10);

  prot::Argument argu_processor;

  bool success = argu_processor.parse(argc, argv);

  if (!success) {
    return 1;
  }

  std::time_t start = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

  std::string start_time_bak = buf;

  std::map<std::string, std::string> arguments = argu_processor.getArguments();

  std::string base_name = arguments["combinedOutputName"];

  std::vector<std::string> spec_file_lst = argu_processor.getSpecFileList();

  std::sort(spec_file_lst.begin(), spec_file_lst.end());

  std::cout << "TopPIC " << prot::version_number << std::endl;

  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    arguments["spectrumFileName"] = spec_file_lst[k];
    prot::TopPICProgress(arguments);
  }

  if (spec_file_lst.size() > 1) {
    arguments["start_time"] = start_time_bak;
    std::cout << "Merging files - started." << std::endl;
    int N = 100000;
    // merge msalign files
    prot::msalign_util::merge_msalign_files(spec_file_lst, N, base_name + "_ms2.msalign");
    // merge feature files
    std::vector<std::string> feature_file_lst(spec_file_lst.size());
    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      std::string sp_file_name = spec_file_lst[i];
      feature_file_lst[i] = sp_file_name.substr(0, sp_file_name.length() - 12) + ".feature";
    }
    prot::feature_util::merge_feature_files(feature_file_lst, N, base_name + ".feature");
    // merge EVALUE files
    std::vector<std::string> prsm_file_lst(spec_file_lst.size());
    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      prsm_file_lst[i] = prot::file_util::basename(spec_file_lst[i]) + ".EVALUE"; 
    }
    prot::prsm_util::merge_prsm_files(prsm_file_lst, N, base_name + "_ms2.EVALUE");
    std::cout << "Merging files - finished." << std::endl;

    std::string sp_file_name = base_name + "_ms2.msalign";
    arguments["spectrumFileName"] = sp_file_name;

    prot::TopPIC_post(arguments);
  }

  if (arguments["keepTempFiles"] != "true") {
    std::cout << "Deleting temporary files - started." << std::endl;
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      std::string sp_file_name = spec_file_lst[k];
      prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_proteoform_cutoff_xml");
      prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_prsm_cutoff_xml");
      prot::file_util::cleanDir(ori_db_file_name, sp_file_name);
    }

    std::string sp_file_name = base_name + "_ms2.msalign";
    prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_proteoform_cutoff_xml");
    prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_prsm_cutoff_xml");
    prot::file_util::cleanDir(ori_db_file_name, sp_file_name);

    std::cout << "Deleting temporary files - finished." << std::endl; 
  }

  return EXIT_SUCCESS;
}
