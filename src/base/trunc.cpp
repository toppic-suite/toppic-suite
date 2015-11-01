#include "base/logger.hpp"
#include "base/acid_util.hpp"
#include "base/trunc.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Trunc::Trunc(const std::string &name, int trunc_len, 
             const std::string &acid_str) {
  name_ = name;
  trunc_len_ = trunc_len;
  acid_ptr_vec_ = AcidUtil::convertStrToAcidPtrVec(acid_str);
  shift_ = -AcidUtil::compAcidPtrVecMass(acid_ptr_vec_);
}

Trunc::Trunc(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  trunc_len_ = XmlDomUtil::getIntChildValue(element, "trunc_len", 0);
  std::string acid_str = XmlDomUtil::getChildValue(element, "acid_str", 0);
  LOG_DEBUG( "name " << name_ << " str " << acid_str << " trunc len " << trunc_len_);
  acid_ptr_vec_ = AcidUtil::convertStrToAcidPtrVec(acid_str);
  shift_ = -AcidUtil::compAcidPtrVecMass(acid_ptr_vec_);
}

std::string Trunc::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

}
