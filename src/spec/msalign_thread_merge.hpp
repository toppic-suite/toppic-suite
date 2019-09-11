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

#ifndef TOPPIC_SPEC_MSALIGN_THREAD_MERGE_HPP_
#define TOPPIC_SPEC_MSALIGN_THREAD_MERGE_HPP_

#include <memory>
#include <vector>
#include <string>

namespace toppic {

class MsalignThreadMerge {
 public:
  MsalignThreadMerge(const std::string &spec_file_name,
                     const std::string &in_file_ext,
                     int in_num,
                     const std::string &out_file_ext);

  void process();

 private:
  std::string spec_file_name_;
  std::vector<std::string> input_file_exts_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<MsalignThreadMerge> MsalignThreadMergePtr;
} /* namespace toppic */

#endif
