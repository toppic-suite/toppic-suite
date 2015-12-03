#ifndef PROT_BASE_PROT_MOD_BASE_HPP_
#define PROT_BASE_PROT_MOD_BASE_HPP_

#include "base/ptm.hpp"
#include "base/trunc.hpp"
#include "base/prot_mod.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class ProtModBase {
 public:
  static void initBase(const std::string &file_name);

  static const ProtModPtrVec& getBaseProtModPtrVec() {
    return prot_mod_ptr_vec_;
  }

  static ProtModPtr getProtModPtrByName (const std::string &name);

  static ProtModPtr getProtModPtr_NONE () {
    return prot_mod_ptr_NONE_;
  }

  static ProtModPtr getProtModPtrFromXml(xercesc::DOMElement * element);

  /*
  static ProtModPtr getProtModPtr_NME () {
    return prot_mod_ptr_NME_;
  }

  static ProtModPtr getProtModPtr_NME_ACETYLATION () {
    return prot_mod_ptr_NME_ACETYLATION_;
  }
  */

 private:
  static ProtModPtrVec prot_mod_ptr_vec_;
  static ProtModPtr prot_mod_ptr_NONE_;
  //static ProtModPtr prot_mod_ptr_NME_;
  //static ProtModPtr prot_mod_ptr_NME_ACETYLATION_;

  static std::string getName_NONE() {return "NONE";}
  //static std::string getName_NME() {return "NME";}
  //static std::string getName_NME_ACETYLATION() {return "NME_ACETYLATION";}
};

}
#endif
