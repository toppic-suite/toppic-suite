#ifndef PROT_BASE_ION_TYPE_BASE_HPP_
#define PROT_BASE_ION_TYPE_BASE_HPP_

#include "base/ion_type.hpp"

namespace prot {

class IonTypeBase {
 public:
  static void initBase(const std::string &file_name);
  static IonTypePtrVec getBaseIonTypePtrVec() {return ion_type_ptr_vec_;}
  static IonTypePtr getIonTypePtrByName(const std::string &name);
  static IonTypePtr getIonTypePtr_PREC() {return ion_type_ptr_PREC_;}
  static IonTypePtr getIonTypePtr_B() {return ion_type_ptr_B_;}

  static std::string getName_B() {return "B";}
  static std::string getName_PREC() {return "PREC";}

 private:
  static IonTypePtrVec ion_type_ptr_vec_;
  static IonTypePtr ion_type_ptr_B_;
  static IonTypePtr ion_type_ptr_PREC_;
};

}

#endif
