#ifndef PROT_BASE_ACTIVATION_BASE_HPP_
#define PROT_BASE_ACTIVATION_BASE_HPP_

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

  static ActivationPtr getActivationPtrFromXml(xercesc::DOMElement * element);

};

} /* namespace prot */

#endif /* ACTIVATION_HPP_ */
