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
      xercesc::DOMElement* root = doc->getDocumentElement();
      xercesc::DOMElement* conf = getChildElement(root, "configuration", 0);
      std::string acid_file_name = getChildValue(conf, "acid_list_file_name", 0);
      acid_list_ = getAcidPtrVecInstance(acid_file_name.c_str());
      std::string ptm_file_name = getChildValue(conf, "ptm_list_file_name", 0);
      ptm_list_ = getPtmPtrVecInstance(ptm_file_name.c_str());
      std::string residue_file_name = getChildValue(conf, "residue_list_file_name", 0);
      residue_list_ = getResiduePtrVecInstance(acid_list_,
                                               ptm_list_,
                                               residue_file_name.c_str());
    }
    delete doc;
  }
  delete parser;
}

}

