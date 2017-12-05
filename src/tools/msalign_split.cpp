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
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>

#include "base/file_util.hpp"
#include "spec/ms_header.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"

using namespace prot;

int main(int argc, const char *argv[]) {
  try {
    namespace po = boost::program_options;
    po::options_description desc("Options");
    int N = 10;
    int group_spec_num = 1;
    std::string msalign_file;
    desc.add_options()
        ("help,h", "Help message")
        ("num,N", po::value<int>(&N), "The number of msalign files to split into. Default value: 10.")
        ("group-num,G", po::value<int>(&group_spec_num), "The group number in the msalign file. Default value: 1.")
        ("msalign,m", po::value<std::string>(&msalign_file)->required(), "The msalign file to split");
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
    if (vm.count("help")) {
      std::cout << desc << "\n";
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    if (vm.count("msalign")) {
      if (!boost::filesystem::exists(msalign_file)) {
        std::cerr << msalign_file << " does not exist!" << std::endl;
        return EXIT_FAILURE;
      }
    }

    MsAlignReader ms_reader(msalign_file, group_spec_num, nullptr, std::set<std::string>());

    std::vector<std::ofstream *> output_vec;
    for (int i = 0; i < N; i++) {
      output_vec.push_back(new std::ofstream(prot::file_util::basename(msalign_file) + "_" + std::to_string(i) + ".msalign"));
    }

    std::vector<std::string> ms_lines = ms_reader.readOneSpectrum();
    int cnt = 0;
    while (ms_lines.size() > 0) {
      cnt = cnt % N;
      for (int i = 0; i < group_spec_num; i++) {
        (*output_vec[cnt]) << std::endl;
        for (size_t k = 0; k < ms_lines.size(); k++) {
          (*output_vec[cnt]) << ms_lines[k] << std::endl;
        }
        (*output_vec[cnt]) << std::endl;
        ms_lines = ms_reader.readOneSpectrum();
      }
      cnt++;
    }

    ms_reader.close();

    for (size_t i = 0; i < output_vec.size(); i++) {
      output_vec[i]->close();
    }

  } catch (const std::exception &ex) {
    std::cerr << ex.what() << std::endl;
  }
  return EXIT_SUCCESS;
}

