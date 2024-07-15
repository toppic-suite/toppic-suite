//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "para/sp_para.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_frac_merge.hpp"

namespace toppic {

namespace msalign_frac_merge {

void mergeMsalignFiles(const std::vector<std::string> &msalign_file_names,
                       const std::string &output_file, 
                       const std::string &para_str) {
  std::ofstream outfile; 
  outfile.open(output_file);
  outfile << para_str;

  for (size_t i = 0; i < msalign_file_names.size(); i++) {
    MsAlignReader sp_reader(msalign_file_names[i]); 
    std::vector<std::string> ms_lines = sp_reader.readOneStrSpectrum();
    while (ms_lines.size() > 0) {
      for (size_t k = 0; k< ms_lines.size(); k++) {
        if (ms_lines[k].substr(0, 12) == "SPECTRUM_ID=") {
          outfile << "SPECTRUM_ID=" << (SpPara::getMaxSpecNumPerFile() * i + std::stoi(ms_lines[k].substr(12))) 
            << std::endl;
        } else if (ms_lines[k].substr(0, 10) == "MS_ONE_ID=") {
          outfile << "MS_ONE_ID=" 
            << (SpPara::getMaxSpecNumPerFile() * i + std::stoi(ms_lines[k].substr(10))) << std::endl;
        } else if (ms_lines[k].substr(0, 21) == "PRECURSOR_FEATURE_ID=") {
          outfile << "PRECURSOR_FEATURE_ID=";
          std::string data = ms_lines[k].substr(21);
          if (data != "") {
            std::vector<std::string> ids = str_util::split(data, ":");
            outfile << (SpPara::getMaxFeatureNumPerFile() * i + std::stoi(ids[0])); 
            for (size_t cnt = 1; cnt < ids.size(); cnt++) {
              outfile << ":" << (SpPara::getMaxFeatureNumPerFile() * i + std::stoi(ids[cnt])); 
            }
          }
          outfile << std::endl;
        } else {
          outfile << ms_lines[k] << std::endl;
        }
      }
      outfile << std::endl;
      ms_lines = sp_reader.readOneStrSpectrum();
    }
  }
  outfile.close();
}

void mergeFractions(const std::vector<std::string> &spec_file_names,
                    const std::string &output_file_name,
                    std::string &para_str) {
  std::vector<std::string> ms1_file_names;
  std::vector<std::string> ms2_file_names;
  for (size_t i = 0; i < spec_file_names.size(); i++) { 
    std::string base_name = file_util::basename(spec_file_names[i]);
    std::string ms1_name = base_name + "_ms1.msalign";
    ms1_file_names.push_back(ms1_name);
    std::string ms2_name = base_name + "_ms2.msalign";
    ms2_file_names.push_back(ms2_name);
  }

  std::string ms1_spec_output_name = output_file_name + "_ms1.msalign";
  std::string ms2_spec_output_name = output_file_name + "_ms2.msalign";

  mergeMsalignFiles(ms1_file_names, ms1_spec_output_name, para_str); 
  mergeMsalignFiles(ms2_file_names, ms2_spec_output_name, para_str); 
}

}

} // namespace toppic 
