#include <iostream>
#include <iomanip>

#include "feature/deconv_argu.hpp"
#include "feature/deconv_para.hpp"
#include "feature/deconv_process.hpp"
#include "feature/feature_detect.hpp"

namespace prot {

int deconvProcess(int argc, char* argv[]) {
  try {
    /*
    DeconvArgument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    LOG_DEBUG("parse complete");
    DeconvParaPtr para_ptr(new DeconvPara(arguments));
    LOG_DEBUG("deconv para");
    DeconvProcess process(para_ptr);
    LOG_DEBUG("init process");
    process.process();
    */
    std::string file_name(argv[1]);
    FeatureDetect::process(file_name);    
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopFD finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  prot::log_level = 2;
  return prot::deconvProcess(argc, argv);
}
