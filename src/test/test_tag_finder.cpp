#include <iostream>
#include <iomanip>

#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_para.hpp"

#include "tagfinder/tagfinder.hpp"

#include "console/argument.hpp"

namespace prot {

int two_base_opt(int argc, char* argv[]) {
  try {
    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopPIC 1.0" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    argu_processor.outputArguments();

    BaseData::init(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string sp_file_name = arguments["spectrumFileName"];

    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));
    LOG_DEBUG("prsm para inited");

    MsAlignUtil::geneSpIndex(sp_file_name);

    std::cout << "Tag searching started." << std::endl;
    TagfinderPtr tagfinder = std::make_shared<Tagfinder>(prsm_para_ptr, 3);
    tagfinder->process();
    tagfinder = nullptr;
    std::cout << "Tag searching finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPC finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  //prot::log_level = 2;
  std::cout << std::setprecision(10);
  return prot::two_base_opt(argc, argv);
}
