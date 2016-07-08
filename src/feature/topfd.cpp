#include <iostream>
#include <iomanip>

#include "feature/deconv_argu.hpp"
#include "feature/deconv_para.hpp"
#include "feature/deconv_process.hpp"

namespace prot {

int deconvProcess(int argc, char* argv[]) {
  try {
    DeconvArgument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    DeconvParaPtr para_ptr(new DeconvPara(arguments));
    DeconvProcess process(para_ptr);
    process.process();
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPC finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  return prot::deconvProcess(argc, argv);
}
