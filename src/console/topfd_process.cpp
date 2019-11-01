//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
//
#include <iomanip>

#include "common/util/logger.hpp"
#include "common/util/time_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/base/base_data.hpp"
#include "spec/msalign_frac_merge.hpp"
#include "feature_detect/feature_detect.hpp"
#include "feature/feature_merge.hpp"
#include "deconv/env/env_base.hpp"
#include "deconv/deconv/deconv_process.hpp"
#include "deconv/deconv/deconv_json_merge.hpp"
#include "console/topfd_argument.hpp"

namespace toppic {

namespace topfd_process {

int processOneFile(std::map<std::string, std::string> arguments, 
                   std::string argu_str, 
                   const std::string &spec_file_name, 
                   int frac_id) {
  try {
    std::cout << "Processing " << spec_file_name << " started." << std::endl;
    std::cout << "Deconvolution started." << std::endl;
    DeconvParaPtr para_ptr = std::make_shared<DeconvPara>(arguments, argu_str);
    DeconvProcess processor(para_ptr, argu_str, spec_file_name, frac_id);
    processor.process();
    std::cout << "Deconvolution finished." << std::endl;

    std::cout << "Feature detection started." << std::endl;
    std::string resource_dir = arguments["resourceDir"];
    bool missing_level_one = (arguments["missingLevelOne"] == "true");
    feature_detect::process(frac_id, spec_file_name, resource_dir,
                            missing_level_one, argu_str);
    std::cout << "Feature detection finished." << std::endl;
    std::cout << "Processing " << spec_file_name << " finished." << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}


void moveFiles(std::string &spec_file_name, bool move_mzrt) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string file_dir =  base_name + "_file";
  file_util::createFolder(file_dir);
  std::string file_name = base_name + "_ms1.msalign";
  file_util::moveFile(file_name, file_dir);
  file_name = base_name + "_frac.feature";
  file_util::moveFile(file_name, file_dir);
  if (move_mzrt) {
    file_name = base_name + "_frac.mzrt.csv";
    file_util::moveFile(file_name, file_dir);
  }
  /*
  if (move_sample_feature) {
    file_name = base_name + "_ms1.feature";
    file_util::moveFile(file_name, file_dir);
  }
  */
}

int process(std::map<std::string, std::string> arguments, 
            std::vector<std::string> spec_file_lst) {
  base_data::init();
  std::string print_str = Argument::geneArgumentStr(arguments, "");
  time_util::addTimeStamp(print_str);
  std::cout << print_str;

  std::string argument_str = Argument::geneArgumentStr(arguments, "#");
  std::string resource_dir = arguments["resourceDir"];
  EnvBase::initBase(resource_dir);
  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    if (str_util::endsWith(spec_file_lst[k], "mzML")
        || str_util::endsWith(spec_file_lst[k], "mzXML")
        || str_util::endsWith(spec_file_lst[k], "mzml")
        || str_util::endsWith(spec_file_lst[k], "mzxml")) {
      int result = processOneFile(arguments, argument_str, spec_file_lst[k], k);
      if (result != 0) {
        return 1;
      }
    }
  }

  if (spec_file_lst.size() > 1) {
    time_util::addTimeStamp(argument_str);
    std::string sample_name = arguments["sampleName"];
    std::cout << "Merging files started." << std::endl;
    MsAlignFracMergePtr msalign_merger = std::make_shared<MsAlignFracMerge>(spec_file_lst, sample_name);
    msalign_merger->process(argument_str);
    msalign_merger = nullptr;
    DeconvJsonMergePtr json_merger = std::make_shared<DeconvJsonMerge>(spec_file_lst, sample_name);
    json_merger->process();
    json_merger = nullptr;
    FeatureMergePtr feature_merger = std::make_shared<FeatureMerge>(spec_file_lst, sample_name);
    feature_merger->process(argument_str);
    feature_merger = nullptr;
    std::cout << "Merging files finished." << std::endl;
  }


  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    if (str_util::endsWith(spec_file_lst[k], "mzML")
        || str_util::endsWith(spec_file_lst[k], "mzXML")
        || str_util::endsWith(spec_file_lst[k], "mzml")
        || str_util::endsWith(spec_file_lst[k], "mzxml")) {
      bool move_mzrt = true;
      moveFiles(spec_file_lst[k], move_mzrt); 
    }
  }

  if (spec_file_lst.size() > 1) {
    std::string sample_name = arguments["sampleName"];
    bool move_mzrt = false;
    moveFiles(sample_name, move_mzrt);
  }

  std::cout << "TopFD finished." << std::endl << std::flush;
  return 0;
}


} // namespace topfd_process 

}  // namespace toppic

