#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/ptm_base.hpp"
#include "local_anno.hpp"

namespace prot {

LocalAnno::LocalAnno(xercesc::DOMElement* element) {
  conf_ = XmlDomUtil::getDoubleChildValue(element, "confidence", 0);
  std::string scr_str = XmlDomUtil::getChildValue(element, "score_list", 0);
  std::vector<std::string> tmp = StringUtil::split(scr_str, ' ');
  for (size_t i = 0; i < tmp.size(); i++) {
    scr_vec_.push_back(std::stod(tmp[i]));
  }
  std::string ptm_element_name = Ptm::getXmlElementName();
  int ptm_count = XmlDomUtil::getChildCount(element, ptm_element_name.c_str());
  if (ptm_count == 0) {
    ptm_ptr_ = nullptr;
  } else {
    xercesc::DOMElement* ptm_element 
        = XmlDomUtil::getChildElement(element, ptm_element_name.c_str(), 0);
    ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);        
  }
}

void LocalAnno::appendToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent) {
  std::string element_name = getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(conf_);
  xml_doc->addElement(element, "confidence", str.c_str());

  str = StringUtil::convertToString(scr_vec_[0]);
  for (size_t i = 1; i < scr_vec_.size(); i++) {
    str = str + " " + StringUtil::convertToString(scr_vec_[i]);
  }

  xml_doc->addElement(element, "score_list", str.c_str());

  if (ptm_ptr_ != nullptr) {
    ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
  }
  parent->appendChild(element);
}

}
