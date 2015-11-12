#ifndef PROT_BASE_MOD_HPP_
#define PROT_BASE_MOD_HPP_

#include "base/residue.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Mod {
 public:
  Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr);

  Mod(xercesc::DOMElement* element); 

  ResiduePtr getOriResiduePtr() { return ori_residue_ptr_;}

  ResiduePtr getModResiduePtr() { return mod_residue_ptr_;}

  double getShift() {return mod_residue_ptr_->getMass() - ori_residue_ptr_->getMass();}

  void appendToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "mod";}

 private:
  ResiduePtr ori_residue_ptr_;
  ResiduePtr mod_residue_ptr_;
};

typedef std::shared_ptr<Mod> ModPtr;
typedef std::vector<ModPtr> ModPtrVec;

}
#endif
