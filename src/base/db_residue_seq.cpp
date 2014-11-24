#include "base/db_residue_seq.hpp"

namespace prot {

DbResidueSeq::DbResidueSeq(const ResiduePtrVec &residues, int id, 
                           const std::string &name,
                           const std::string &desc)
    : ResidueSeq(residues) {
      id_ = id;
      name_ = name;
      desc_ = desc;
    }

void DbResidueSeq::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("db_residue_seq");
  std::string str = convertToString(id_);
  xml_doc->addElement(element, "db_seq_id", str.c_str());
  xml_doc->addElement(element, "db_seq_name", name_.c_str());
  xml_doc->addElement(element, "db_seq_desc", desc_.c_str());
  parent->appendChild(element);
}

}
