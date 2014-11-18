#ifndef PROT_DB_RESIDUE_SEQ_HPP_
#define PROT_DB_RESIDUE_SEQ_HPP_

#include <string>

#include "base/residue_seq.hpp"

namespace prot {

/* database residue sequence contains protein name and id */
class DbResidueSeq: public ResidueSeq {
 public:
  DbResidueSeq(const ResiduePtrVec &residues, int id, 
               const std::string &name,
               const std::string &desc);

  int getId() {return id_;}
  const std::string& getName() {return name_;}
  const std::string& getDesc() {return desc_;}
  std::string getNameDesc() {return name_ + " " + desc_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  int id_;
  std::string name_;
  std::string desc_;
};

typedef std::shared_ptr<DbResidueSeq> DbResSeqPtr;

}
#endif
