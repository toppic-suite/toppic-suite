//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>

#include "spec/ms_header.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/msalign_reader.hpp"

int main(int argc, const char *argv[]) {
  try {
    namespace po = boost::program_options;
    po::options_description desc("Options");
    int N = 100000;
    std::string output_filename = "";
    std::vector<std::string> msalign_file_list;
    desc.add_options()
        ("help,h", "Help message")
        ("num,N", po::value<int>(&N), "N")
        ("msalign,m", po::value<std::vector<std::string> >()->multitoken()->required(), "msalign file list")
        ("output,o", po::value<std::string>(&output_filename)->required(), "output msalign file name");
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
    if (vm.count("help")) {
      std::cout << desc << "\n";
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    if (vm.count("msalign")) {
      msalign_file_list = vm["msalign"].as<std::vector<std::string> >();
      if (msalign_file_list.size() == 1) {
        return EXIT_SUCCESS;
      }
    }

    std::cout << "msalign files to merge:" << std::endl;
    for (size_t i = 0; i < msalign_file_list.size(); i++) {
      if (!boost::filesystem::exists(msalign_file_list[i])) {
        std::cerr << "msalign file " << msalign_file_list[i] << " does not exist!" << std::endl;
        return EXIT_FAILURE;
      }
      std::cout << msalign_file_list[i] << std::endl;
    }

    std::ofstream outfile(output_filename); 

    for (int i = 0; i < static_cast<int>(msalign_file_list.size()); i++) {
      prot::MsAlignReader sp_reader(msalign_file_list[i], 1, nullptr, std::set<std::string>());
      std::vector<std::string> ms_lines = sp_reader.readOneSpectrum();
      while (ms_lines.size() > 0) {
        for (size_t k = 0; k < ms_lines.size(); k++) {
          if (ms_lines[k].substr(0, 3) == "ID=") {
            outfile << "ID=" << (N * i + std::stoi(ms_lines[k].substr(3))) << std::endl;
          } else if (ms_lines[k].substr(0, 10) == "MS_ONE_ID=") {
            outfile << "MS_ONE_ID=" << (N * i + std::stoi(ms_lines[k].substr(10))) << std::endl;
          } else if (ms_lines[k].substr(0, 6) == "SCANS=") {
            outfile << "SCANS=" << (N * i + std::stoi(ms_lines[k].substr(6))) << std::endl;
          } else {
            outfile << ms_lines[k] << std::endl;
          }
        }
        ms_lines = sp_reader.readOneSpectrum();
      }
      sp_reader.close();
    }

    outfile.close();

  } catch (const std::exception &ex) {
    std::cerr << ex.what() << std::endl;
  }
  return EXIT_SUCCESS;
}

