#include <base/logger.hpp>

#include "base/trunc.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

TruncPtrVec TruncFactory::trunc_ptr_vec_;

Trunc::Trunc(const std::string &name, int trunc_len, 
             const std::string &str) {
  name_ = name;
  trunc_len_ = trunc_len;
  shift_ = 0;
  for (size_t i = 0; i < str.length(); i++) {
    std::string letter = str.substr(i, 1);
    AcidPtr acid_ptr = AcidFactory::getBaseAcidPtrByOneLetter(letter);
    acid_str_.push_back(acid_ptr);
    shift_ = shift_ - acid_ptr->getMonoMass();
  }
}

void Trunc::appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("truncation");
  xml_doc->addElement(element, "name", name_.c_str());
  std::string str = convertToString(trunc_len_);
  xml_doc->addElement(element, "trunc_len", str.c_str());
  str = convertToString(shift_);
  xml_doc->addElement(element, "shift", str.c_str());
  xercesc::DOMElement* acid_list = xml_doc->createElement("amino_acid_list");
  for(size_t i=0;i<acid_str_.size();i++){
    acid_str_[i]->appendxml(xml_doc,acid_list);
  }
  element->appendChild(acid_list);
  parent->appendChild(element);
}

bool Trunc::isSameTrunc(int len, const ResiduePtrVec& res_ptr_vec) {
    if(trunc_len_ != len){
        return false;
    }
    for(int i=0;i<trunc_len_;i++){
        if(acid_str_[i] != res_ptr_vec[i]->getAcidPtr()){
            return false;
        }
    }
    return true;
}

bool Trunc::isValidTrunc(const ResiduePtrVec & res_ptr_vec) {
    //check if trunc acids match N-terminal acids of the protein 
    bool result = true;
    if (trunc_len_ >= (int)res_ptr_vec.size()) {
        result = false; 
    }
    else {
        result = isSameTrunc(trunc_len_, res_ptr_vec);
    }
    //LOG_DEBUG("Valid trunc " << result << " trunc len " << trunc_len_ 
    //          << " seq len " << res_seq_ptr->getLen());
    return result;
}

void TruncFactory::initFactory(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int trunc_num = getChildCount(parent, "truncation");
    for (int i = 0; i < trunc_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "truncation", i);
      std::string name = getChildValue(element, "name", 0);
      int trunc_len = getIntChildValue(element, "trunc_len", 0);
      std::string str = getChildValue(element, "acid_str", 0);
      LOG_DEBUG( "name " << name << " str " << str << " trunc len " << trunc_len);
      trunc_ptr_vec_.push_back(TruncPtr(new Trunc(name, trunc_len, str)));
    }
  }
}

TruncPtr TruncFactory::getBaseTruncPtrByName(const std::string &name) {
  for (size_t i = 0; i < trunc_ptr_vec_.size(); i++) {
    std::string n = trunc_ptr_vec_[i]->getName();
    if (n == name) {
      return trunc_ptr_vec_[i];
    }
  }
  return TruncPtr(nullptr);
}

}
