#ifndef PROT_BASE_TRUNC_BASE_HPP_
#define PROT_BASE_TRUNC_BASE_HPP_

#include <string>
#include "base/trunc.hpp"

namespace prot {

class TruncBase {
 public:
  static void initBase(const std::string &file_name);
  static const TruncPtrVec& getBaseTruncPtrVec() {return trunc_ptr_vec_;}
  static TruncPtr getTruncPtrByName(const std::string &name);

 private:
  static TruncPtrVec trunc_ptr_vec_;
};

}

#endif
