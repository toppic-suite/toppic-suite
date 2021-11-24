//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "topfd/common/topfd_para.hpp"
#include "topfd/common/topfd_process.hpp"
#include "console/topfd_argument.hpp"

int main(int argc, char* argv[]) {
  toppic::Argument argu_processor;
  bool success = argu_processor.parse(argc, argv);

  if (!success) {
    return 1;
  }

  toppic::TopfdParaPtr topfd_para_ptr = argu_processor.getTopfdParaPtr();

  std::vector<std::string> spec_file_lst = argu_processor.getSpecFileList();

  int result = toppic::topfd_process::process(topfd_para_ptr, spec_file_lst);

  return result;
}

