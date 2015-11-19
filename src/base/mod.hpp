#ifndef PROT_BASE_MOD_HPP_
#define PROT_BASE_MOD_HPP_

#include "base/residue.hpp"
#include "base/residue_base.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Mod;
typedef std::shared_ptr<Mod> ModPtr;

class Mod {
 public:
  Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr);

  Mod(xercesc::DOMElement* element); 

  ResiduePtr getOriResiduePtr() { return ori_residue_ptr_;}

  ResiduePtr getModResiduePtr() { return mod_residue_ptr_;}

  bool isSame(ModPtr mod_ptr) {
    return ori_residue_ptr_ == mod_ptr->getOriResiduePtr() 
        && mod_residue_ptr_ == mod_ptr->getModResiduePtr();
  }

  double getShift() {return mod_residue_ptr_->getMass() - ori_residue_ptr_->getMass();}

  void appendToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "mod";}

 private:
  ResiduePtr ori_residue_ptr_;
  ResiduePtr mod_residue_ptr_;
};

typedef std::vector<ModPtr> ModPtrVec;

}
#endif
