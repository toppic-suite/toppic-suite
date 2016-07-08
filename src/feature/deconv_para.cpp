#include "feature/deconv_para.hpp"

namespace prot {

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
