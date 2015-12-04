#include "base/logger.hpp"
#include "base/trunc.hpp"
#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Trunc::Trunc(const std::string &name, int trunc_len, 
             const std::string &residue_str) {
  name_ = name;
  trunc_len_ = trunc_len;
  residue_ptr_vec_ = ResidueUtil::convertStrToResiduePtrVec(residue_str);
  shift_ = -ResidueUtil::compResiduePtrVecMass(residue_ptr_vec_);
}

Trunc::Trunc(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  trunc_len_ = XmlDomUtil::getIntChildValue(element, "trunc_len", 0);
  std::string residue_str = XmlDomUtil::getChildValue(element, "residue_str", 0);
  LOG_DEBUG( "name " << name_ << " str " << residue_str << " trunc len " << trunc_len_);
  residue_ptr_vec_ = ResidueUtil::convertStrToResiduePtrVec(residue_str);
  shift_ = -ResidueUtil::compResiduePtrVecMass(residue_ptr_vec_);
}

std::string Trunc::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

}
