#ifndef PROT_BASE_MOD_READER_HPP_
#define PROT_BASE_MOD_READER_HPP_

#include "base/mod.hpp"

namespace prot {

class ModReader {
 public:
  static ModPtrVec readMod(const std::string &file_name);
};

}
#endif
