/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_RESIDUE_BASE_HPP_
#define PROT_BASE_RESIDUE_BASE_HPP_

#include <string>
#include <memory>
#include <map>

#include "base/residue.hpp"
#include "base/logger.hpp"

namespace prot {

class ResidueBase {
 public:
  static void initBase(const std::string &file_name);

  static const ResiduePtrVec& getBaseResiduePtrVec() {return residue_ptr_vec_;}

  static ResiduePtr getEmptyResiduePtr() {return empty_residue_ptr_;}

  static ResiduePtr getResiduePtrFromXml(xercesc::DOMElement * element);
  
  //static ResiduePtrVec getResiduePtrVecInstance(const std::string &file_name);
  
 private:
  static ResiduePtrVec residue_ptr_vec_;
  static ResiduePtr empty_residue_ptr_;

  static ResiduePtr getBaseResiduePtr(ResiduePtr residue_ptr);
};

}

#endif
