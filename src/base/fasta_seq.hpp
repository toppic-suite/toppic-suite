#ifndef PROT_BASE_FASTA_SEQ_HPP_
#define PROT_BASE_FASTA_SEQ_HPP_

#include <memory>
#include <string>
#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class FastaSeq {
 public:
  FastaSeq(const std::string &name_line, const std::string &ori_seq);

  FastaSeq(const std::string &name, const std::string &desc, 
           const std::string &ori_seq);

  std::string getName() {return name_;}

  std::string getDesc() {return desc_;}

  std::string getSeq() {return seq_;}

  int getLen() {return seq_.length();}

  static std::string getXmlElementName() {return "fasta_seq";}

  void appendNameDescToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getNameFromXml(xercesc::DOMElement* element);

  static std::string getDescFromXml(xercesc::DOMElement* element);

 private:
  std::string name_;
  std::string desc_;
  std::string seq_;

  static std::string rmChar(const std::string &ori_seq);
}; 

typedef std::shared_ptr<FastaSeq> FastaSeqPtr;

}

#endif
