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

#include "ms/spec/msalign_reader_util.hpp"

namespace toppic {

namespace msalign_reader_util {

void readAllSpectra(const std::string &msalign_file_name, 
                    DeconvMsPtrVec &ms_ptr_vec) {
  MsAlignReader sp_reader(msalign_file_name); 

  DeconvMsPtr ms_ptr;
  while ((ms_ptr = sp_reader.getNextMsPtr())!= nullptr) {
    ms_ptr_vec.push_back(ms_ptr);
  }
}

}

}  // namespace toppic
