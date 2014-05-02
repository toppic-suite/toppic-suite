#include <iostream>

#include "base/base_data.hpp"
#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

int main(int argc, char* argv[]) {
  try {
    std::cout << "EValueConsole 0.1 " << std::endl;
    prot::TdgfMngPtr mng_ptr = 
        prot::TdgfMngPtr(new prot::TdgfMng (argv[1], argv[2], argv[3], argv[4], argv[5]));
    prot::EValueProcessor processor(mng_ptr);
    processor.init();
    /* compute E-value for a set of prsms each run */
    processor.process(false);
  } catch (const char* e) {
    std::cout << "Exception " << e << std::endl;
  }
  return 0;
}
