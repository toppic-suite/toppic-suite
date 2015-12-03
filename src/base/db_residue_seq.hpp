#ifndef PROT_BASE_DB_RESIDUE_SEQ_HPP_
#define PROT_BASE_DB_RESIDUE_SEQ_HPP_

#include <string>

#include "base/residue_seq.hpp"

namespace prot {

/* database residue sequence contains protein name and id */
class DbResidueSeq: public ResidueSeq {
 public:
  DbResidueSeq(const ResiduePtrVec &residues,  
               const std::string &name,
               const std::string &desc);

  const std::string& getName() {return name_;}
  const std::string& getDesc() {return desc_;}
  std::string getNameDesc() {return name_ + " " + desc_;}

  static std::string getXmlElementName() {return "db_residue_seq";}

  void appendRefToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  std::string name_;
  std::string desc_;
};

typedef std::shared_ptr<DbResidueSeq> DbResSeqPtr;

}
#endif
