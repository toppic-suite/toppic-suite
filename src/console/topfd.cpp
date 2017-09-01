// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <iostream>
#include <iomanip>
#include <map>
#include <string>

#include "base/file_util.hpp"
#include "base/base_data.hpp"
#include "feature/deconv_argu.hpp"
#include "feature/deconv_para.hpp"
#include "feature/deconv_process.hpp"
#include "feature/feature_detect.hpp"

namespace prot {

int deconvProcess(int argc, char* argv[]) {
  try {
    time_t start = time(0);
    std::string exe_dir = FileUtil::getExecutiveDir(argv[0]);
    BaseData::init(exe_dir);
    DeconvArgument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    LOG_DEBUG("parse complete");

    DeconvParaPtr para_ptr = std::make_shared<DeconvPara>(arguments);
    LOG_DEBUG("deconv para");
    DeconvProcess process(para_ptr);
    LOG_DEBUG("init process");
    process.process();
    FeatureDetect::process(para_ptr);

    time_t end = time(0);
    std::cout << "Runing time: " << std::to_string(static_cast<int>(difftime(end, start))) << " seconds." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << std::endl << "TopFD finished." << std::endl;
  return 0;
}

}  // namespace prot

int main(int argc, char* argv[]) {
  // prot::log_level = 2;
  return prot::deconvProcess(argc, argv);
}
