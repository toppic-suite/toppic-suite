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

//#include <algorithm>

#include "console/topindex_process.hpp"
#include "gui/topindex/threadtopindex.h"
#include "gui/util/run_exe.h"

void ThreadTopIndex::run() {
  //std::sort(spec_file_lst_.begin(), spec_file_lst_.end());
  //toppic::TopIndexProcess(arguments_);
  toppic::RunExe runExe;
  std::string cmd = runExe.geneCommand(arguments_, "topindex");
  runExe.run(cmd);
}
