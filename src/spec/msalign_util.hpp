#ifndef PROT_SPEC_MSALIGN_UTIL_HPP_
#define PROT_SPEC_MSALIGN_UTIL_HPP_

#include "spec/msalign_reader.hpp"

namespace prot {

class MsAlignUtil {
 public:
  static int countSpNum(const std::string &spectrum_file);

  static void geneSpIndex(const std::string &spectrum_file_name);

  static int getSpNum(const std::string &spectrum_file_name);
};

}
#endif
