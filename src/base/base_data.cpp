/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */

#include "base_data.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

BaseData::BaseData  (const char* config_file_name) {

  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, config_file_name);
    if (doc) {
      xercesc::DOMElement* element = doc->getElement("configuration", 0);
      std::string acid_file_name = getChildValue(element, "acid_list_file_name");
      acid_list_ = getAcidPtrVecInstance(acid_file_name.c_str());
      std::string ptm_file_name = getChildValue(element, "ptm_list_file_name");
      ptm_list_ = getPtmPtrVecInstance(ptm_file_name.c_str());
      std::string residue_file_name = getChildValue(element, "residue_list_file_name");
      residue_list_ = getResiduePtrVecInstance(acid_list_,
                                               ptm_list_,
                                               residue_file_name.c_str());
    }
    delete doc;
  }
  delete parser;
}

}

