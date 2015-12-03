#include "base/logger.hpp"
#include "base/fasta_seq.hpp"

namespace prot {

FastaSeq::FastaSeq(const std::string &name_line, 
                   const std::string &ori_seq) {
  int space_pos = name_line.find(" ");
  name_ = name_line.substr(0, space_pos);
  desc_ = name_line.substr(space_pos + 1);
  seq_ = rmChar(ori_seq);
}

FastaSeq::FastaSeq(const std::string &name, 
                   const std::string &desc, 
                   const std::string &ori_seq): 
    name_(name),
    desc_(desc) {
      seq_ = rmChar(ori_seq);
    }

/** process fasta string and remove unknown letters */
std::string FastaSeq::rmChar(const std::string &ori_seq) {
  std::string seq = "";
  for (size_t i = 0; i < ori_seq.length(); i++) {
    char c = ori_seq.at(i);
    if (c < 'A' || c > 'Z') {
      continue;
    }
    char r = c;
    if (c == 'B') {
      r = 'D';
    } else if (c == 'Z') {
      r = 'E';
    } else if (c == 'X') {
      r = 'A';
    } else if (c == 'J') {
      r = 'I';
    }
    seq = seq + r;
  }
  if (ori_seq != seq) { 
    LOG_INFO( "Reading sequence. Unknown letter occurred. ");
  }
  return seq;
}

void FastaSeq::appendNameDescToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = FastaSeq::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "seq_name", name_.c_str());
  xml_doc->addElement(element, "seq_desc", desc_.c_str());
  parent->appendChild(element);
}


}
