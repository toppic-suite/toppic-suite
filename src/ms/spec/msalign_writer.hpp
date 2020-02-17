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

#ifndef PROT_SPEC_MSALIGN_WRITER_HPP_
#define PROT_SPEC_MSALIGN_WRITER_HPP_

#include <fstream>

#include "ms/spec/deconv_ms.hpp"

namespace toppic {

class MsAlignWriter {
 public:
  MsAlignWriter(const std::string &file_name);

  ~MsAlignWriter();

  void write(DeconvMsPtr ms_ptr);

  void writePara(const std::string &para_str);

  void close();

 private:
  std::string file_name_;
  std::ofstream output_;
};

typedef std::shared_ptr<MsAlignWriter> MsAlignWriterPtr;
typedef std::vector<MsAlignWriterPtr>  MsAlignWriterPtrVec;

}  // namespace toppic

#endif
