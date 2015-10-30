#ifndef PROT_ACTIVATION_BASE_HPP_
#define PROT_ACTIVATION_BASE_HPP_

#include "base/activation.hpp"

namespace prot {

/* activation base */
class ActivationBase {
 private:
  static ActivationPtrVec activation_ptr_vec_;

 public:
  static void initBase(const std::string &file_name);

  static const ActivationPtrVec& getActivationPtrVec() {return activation_ptr_vec_;}

  static ActivationPtr getActivationPtrByName(const std::string &name);
};

} /* namespace prot */

#endif /* ACTIVATION_HPP_ */
