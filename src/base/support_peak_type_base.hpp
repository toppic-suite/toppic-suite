#ifndef PROT_BASE_SUPPORT_PEAK_TYPE_BASE_HPP_
#define PROT_BASE_SUPPORT_PEAK_TYPE_BASE_HPP_

#include <memory>
#include <string>
#include <vector>

#include "base/support_peak_type.hpp"

namespace prot {

class SPTypeBase {
 public:
  static void initBase(const std::string &file_name);

  static const SPTypePtrVec& getBaseSPTypePtrVec() {
    return sp_type_ptr_vec_;}

  static SPTypePtr getSPTypePtrByName(const std::string &name);
  static SPTypePtr getSPTypePtrById(int id);

  static SPTypePtr getSPTypePtr_N_TERM() {
    return sp_type_ptr_N_TERM_;
  }

  static std::string getName_N_TERM() {return "N_TERM";}

 private:
  static SPTypePtrVec sp_type_ptr_vec_;
  static SPTypePtr sp_type_ptr_N_TERM_;
};

} /* namespace prot */

#endif /* SUPPORT_PEAK_TYPE_HPP_ */
