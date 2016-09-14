// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include "base/file_util.hpp"
#include "console/summary.hpp"

namespace prot {

void Summary::outputSummary(Argument arguments, 
                            boost::posix_time::ptime start_time,
                            boost::posix_time::ptime stop_time) {
  
  std::string spec_file_name = arguments.getArguments()["spectrumFileName"];
  std::string base_name = FileUtil::basename(spec_file_name);
  std::string output_file_name = base_name + "." + "SUMMARY";

  std::ofstream output; 
  output.open(output_file_name, std::ios::out | std::ios::app);

  boost::posix_time::time_duration duration;
  duration = stop_time - start_time;
  output << "1. Time" << std::endl;
  output << "Start time: " << boost::posix_time::to_simple_string(start_time) << std::endl;
  output << "Stop time: " << boost::posix_time::to_simple_string(stop_time) << std::endl;
  output << "Running time: " << boost::posix_time::to_simple_string(duration) << std::endl;


  Argument::outputArguments(output, arguments.getArguments());
  output.close();
}

} /* namespace prot */
