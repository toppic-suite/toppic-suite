/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <iostream>

#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

Ptm::Ptm(const std::string &abbr_name, 
            double mono_mass) {
    abbr_name_ = abbr_name;
    mono_mass_ = mono_mass;
}

bool Ptm::isEmpty() {
  if (mono_mass_ == 0.0) {
    return true;
  }
  else {
    return false;
  }
}

PtmPtrVec getPtmPtrVecInstance(const char* file_name) {
  PtmPtrVec ptm_list;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name);
    if (doc) {
      xercesc::DOMElement* parent = doc->getDocumentElement();
      int ptm_num = getChildCount(parent, "modification");
      for (int i = 0; i < ptm_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "modification", i);
        std::string abbr_name = getChildValue(element, "abbreviation", 0);
        double mono_mass = getDoubleChildValue(element, "mono_mass", 0);
        ptm_list.push_back(PtmPtr(new Ptm(abbr_name, mono_mass)));

      }
      delete doc;
    }
  }
  return ptm_list;
}

/**
 *   Returns a PTM based on the abbreviation name. Returns null if the
 *   abbreviation name does not exist.
 */
PtmPtr getPtmPtrByAbbrName(PtmPtrVec &ptm_list, 
                           const std::string &abbr_name) {
  for (unsigned int i = 0; i < ptm_list.size(); i++) {
    std::string n = ptm_list[i]->getAbbrName();
    if (n == abbr_name) {
      return ptm_list[i];
    }
  }
  return PtmPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containAbbrsName(PtmPtrVec &ptm_list, 
                      const std::string &abbr_name) {
    if (getPtmPtrByAbbrName(ptm_list, abbr_name).get() == nullptr) {
        return false;
    }
    else {
        return true;
    }
}

PtmPtr findEmptyPtmPtr(PtmPtrVec &ptm_list) {
  for (unsigned int i = 0; i < ptm_list.size(); i++) {
    if (ptm_list[i]->isEmpty()) {
      return ptm_list[i];
    }
  }
  throw "Empty ptm does not exist!";
}

PtmPtr addPtm(PtmPtrVec &ptm_list, std::string abbr_name,
              double mono_mass) {
  PtmPtr ptm_ptr = getPtmPtrByAbbrName(ptm_list, abbr_name);
  if (ptm_ptr.get() == nullptr) {
    PtmPtr new_ptm(new Ptm(abbr_name, mono_mass));
    ptm_list.push_back(new_ptm);
    return new_ptm;
  }
  else {
    return ptm_ptr;
  }
}

}

