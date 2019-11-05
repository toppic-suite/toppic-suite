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

#include <map>
#include <string>

#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/time_util.hpp"
#include "spec/msalign_frac_merge.hpp"
#include "feature/feature_merge.hpp"
#include "console/topfd_argument.hpp"
#include "console/topfd_process.hpp"

int main(int argc, char* argv[]) {
  toppic::Argument argu_processor;
  bool success = argu_processor.parse(argc, argv);

  if (!success) {
    return 1;
  }

  std::map<std::string, std::string> arguments = argu_processor.getArguments();

  std::string exe_dir = toppic::file_util::getExecutiveDir(argv[0]);

  arguments["executiveDir"] = exe_dir;

  std::vector<std::string> spec_file_lst = argu_processor.getSpecFileList();

  int result = toppic::topfd_process::process(arguments, spec_file_lst);

  return result;
}
