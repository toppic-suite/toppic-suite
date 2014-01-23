#include "base/db_residue_seq.hpp"

namespace prot {

DbResidueSeq::DbResidueSeq(ResiduePtrVec residues, int id, std::string name)
    : ResidueSeq(residues) {
      id_ = id;
      name_ = name;
    }

void DbResidueSeq::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("db_residue_seq");
  std::string str = convertToString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

}
