#include "base/logger.hpp"
#include "base/trunc.hpp"
#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Trunc::Trunc(const std::string &name, int trunc_len, 
             const std::string &trunc_residues,
             const std::string &allow_first_remain_residues) {
  name_ = name;
  trunc_len_ = trunc_len;
  trunc_residue_ptr_vec_ = ResidueUtil::convertStrToResiduePtrVec(trunc_residues);
  allow_first_remain_residue_ptrs_ = ResidueUtil::convertStrToResiduePtrVec(allow_first_remain_residues);
  shift_ = -ResidueUtil::compResiduePtrVecMass(trunc_residue_ptr_vec_);
}

Trunc::Trunc(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  trunc_len_ = XmlDomUtil::getIntChildValue(element, "trunc_len", 0);
  std::string trunc_residues = XmlDomUtil::getChildValue(element, "trunc_residues", 0);
  LOG_DEBUG( "name " << name_ << " str " << trunc_residues << " trunc len " << trunc_len_);
  trunc_residue_ptr_vec_ = ResidueUtil::convertStrToResiduePtrVec(trunc_residues);
  std::string allow_first_remain_residues = XmlDomUtil::getChildValue(element, "allow_first_remain_residues", 0);
  allow_first_remain_residue_ptrs_ = ResidueUtil::convertStrToResiduePtrVec(allow_first_remain_residues);
  shift_ = -ResidueUtil::compResiduePtrVecMass(trunc_residue_ptr_vec_);
}

std::string Trunc::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

}
