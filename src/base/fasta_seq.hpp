#ifndef PROT_BASE_FASTA_SEQ_HPP_
#define PROT_BASE_FASTA_SEQ_HPP_

#include <memory>
#include <string>
#include <vector>

#include "base/residue.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {


class FastaSeq {
 public:
  FastaSeq(const std::string &name_line, const std::string &ori_seq);

  FastaSeq(const std::string &name, const std::string &desc, 
           const std::string &ori_seq);

  std::string getName() {return name_;}

  std::string getDesc() {return desc_;}

  std::string getRawSeq() {return seq_;}

  StringPairVec getAcidPtmPairVec() {return acid_ptm_pair_vec_;}

  int getAcidPtmPairLen() {return acid_ptm_pair_vec_.size();}

  //int getLen() {return seq_.length();}

  static std::string getXmlElementName() {return "fasta_seq";}

  void appendNameDescToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getNameFromXml(xercesc::DOMElement* element);

  static std::string getDescFromXml(xercesc::DOMElement* element);

  static std::string getString(const std::pair<std::string,std::string> &str_pair);

  static std::string getString(const StringPairVec &str_pair_vec);

 private:
  std::string name_;
  std::string desc_;
  std::string seq_;
  StringPairVec acid_ptm_pair_vec_;

  static std::string rmChar(const std::string &ori_seq);
  void compAcidPtmPairVec();
}; 

typedef std::shared_ptr<FastaSeq> FastaSeqPtr;

}

#endif
