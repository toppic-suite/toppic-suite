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

#include <set>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_writer.hpp"
#include "spec/msalign_frac_merge.hpp"

namespace toppic {

namespace msalign_frac_merge {

void mergeFiles(const std::vector<std::string> &spec_file_lst,
                const std::string &output_file, 
                int max_num_per_file,
                const std::string &para_str) {
  std::ofstream outfile; 
  outfile.open(output_file);
  outfile << para_str;

  for (size_t i = 0; i < spec_file_lst.size(); i++) {
    MsAlignReader sp_reader(spec_file_lst[i], 1, nullptr, std::set<std::string>());
    std::vector<std::string> ms_lines = sp_reader.readOneStrSpectrum();
    while (ms_lines.size() > 0) {
      for (size_t k = 0; k< ms_lines.size(); k++) {
        if (ms_lines[k].substr(0, 3) == "ID=") {
          outfile << "ID=" << (max_num_per_file * i + std::stoi(ms_lines[k].substr(3))) 
              << std::endl;
        } else if (ms_lines[k].substr(0, 10) == "MS_ONE_ID=") {
          outfile << "MS_ONE_ID=" 
              << (max_num_per_file * i + std::stoi(ms_lines[k].substr(10))) << std::endl;
        } else {
          outfile << ms_lines[k] << std::endl;
        }
      }
      outfile << std::endl;
      ms_lines = sp_reader.readOneStrSpectrum();
    }
    sp_reader.close();
  }

  outfile.close();
}

}

} /* namespace toppic */
