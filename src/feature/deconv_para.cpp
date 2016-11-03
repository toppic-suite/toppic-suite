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



#include <iostream>
#include "feature/deconv_para.hpp"

namespace prot {

DeconvPara::DeconvPara(std::map<std::string, std::string> &arguments) {
  data_file_name_ = arguments["spectrumFileName"];
  setInputType(arguments["inputType"]);
  setOutputType(arguments["outputType"]);

  ms_level_ = std::stoi(arguments["msLevel"]);
  missing_level_one_ = (arguments["missingLevelOne"] == "true");
  max_charge_ = std::stoi(arguments["maxCharge"]);
  max_mass_ = std::stod(arguments["maxMass"]);
  tolerance_ = std::stod(arguments["mzError"]);
  sn_ratio_ = std::stod(arguments["snRatio"]);
  keep_unused_peaks_ = (arguments["keepUnusedPeaks"] == "true");
  prec_window_ = std::stod(arguments["precWindow"]);
  exec_dir_ = arguments["executiveDir"];
}

int DeconvPara::setOutputType (std::string &format) {
  if (format == "mgf") {
    output_type_ = OUTPUT_MGF;
  } else if (format == "msalign") {
    output_type_ = OUTPUT_MSALIGN;
  } else if (format == "text") {
    output_type_ = OUTPUT_TEXT;
  } else {
    return 1;
  }
  return 0;
}

int DeconvPara::setInputType (std::string &format) {
  if (format=="mgf") {
    input_type_ = INPUT_MGF;
  } else if (format=="mzXML") {
    input_type_ = INPUT_MZXML;
  } else {
    return 1;
  }
  return 0;
}


}
