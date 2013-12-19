#include <iostream>

#include "base/base_data.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

int main(int argc, char* argv[]) {
  try {
    std::cout << "ZeroPtmConsole 0.1 " << std::endl;
    prot::ZeroPtmMngPtr mng_ptr = prot::ZeroPtmMngPtr(new prot::ZeroPtmMng (std::string(argv[1])));
    mng_ptr->search_db_file_name_ = argv[2];
    mng_ptr->spectrum_file_name_ = argv[3];
    mng_ptr->output_file_ext_ = argv[4];
    std::cout << "protein database " << argv[2] << " spectrum dataset " << argv[3] << std::endl;
    prot::zeroPtmSearchProcess(mng_ptr);
    //mng.resFileName = args[1];
    //mng.spPara.setActivationType(EnumActivation.getActivationType(args[4]));
    //ZeroPtmProcessor processor = new ZeroPtmProcessor(mng);
    //processor.process();
  } catch (const char* e) {
    std::cout << "Exception " << e << std::endl;
  } 

  return 0;
}

