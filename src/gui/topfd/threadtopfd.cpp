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

#include <algorithm>

#include "threadtopfd.h"

void ThreadTopFD::run() {
  std::sort(spec_file_lst_.begin(), spec_file_lst_.end());
  for (size_t k = 0; k < spec_file_lst_.size(); k++) {
    if (toppic::str_util::endsWith(spec_file_lst_[k], "mzML")
        || toppic::str_util::endsWith(spec_file_lst_[k], "mzXML")
        || toppic::str_util::endsWith(spec_file_lst_[k], "mzml")
        || toppic::str_util::endsWith(spec_file_lst_[k], "mzxml")) {
      arguments_["spectrumFileName"] = spec_file_lst_[k];
      toppic::TopFDProcess(arguments_);
    }
  }
}
