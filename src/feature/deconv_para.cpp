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
}

int DeconvPara::setOutputType (std::string &format) {
  if (format=="mgf") {
    output_type_ = OUTPUT_MGF;
  } else if (format==("msalign")) {
    output_type_ = OUTPUT_MSALIGN;
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
