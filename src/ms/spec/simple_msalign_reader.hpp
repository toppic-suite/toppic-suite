//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_MS_SPEC_SIMPLE_MSALIGN_READER_HPP_
#define TOPPIC_MS_SPEC_SIMPLE_MSALIGN_READER_HPP_

#include <limits>
#include <fstream>
#include <string>
#include <vector>

#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/spectrum_set.hpp"

namespace toppic {

class SimpleMsAlignReader {
 public:
  SimpleMsAlignReader(const std::string &file_name);

  SimpleMsAlignReader(const std::string &file_name, 
                      int group_spec_num,
                      ActivationPtr activation_ptr);

  ~SimpleMsAlignReader();

  std::vector<std::string> readOneStrSpectrum();

  DeconvMsPtr getNextMsPtr();

  DeconvMsPtrVec getNextMsPtrVec(); 

  static void readMsOneSpectra(const std::string &file_name, DeconvMsPtrVec &ms_ptr_vec); 

 private:
  std::string file_name_;

  int group_spec_num_ = 1;

  ActivationPtr activation_ptr_;

  std::ifstream input_;

  std::vector<std::string> spectrum_str_vec_;

  int current_ = 0;

  DeconvMsPtr deconv_ms_ptr_ = nullptr;

  void readNext();
};

typedef std::shared_ptr<SimpleMsAlignReader> SimpleMsAlignReaderPtr;
typedef std::vector<SimpleMsAlignReaderPtr>  SimpleMsAlignReaderPtrVec;

}  // namespace toppic
#endif
