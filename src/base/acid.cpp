/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <iostream>

#include "acid.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace proteomics {

Acid::Acid (std::string const &name, std::string const &one_letter, 
            std::string const &three_letter, std::string const &composition, 
            double mono_mass, double avg_mass) {
  name_ = name;
  one_letter_ = one_letter;
  three_letter_ = three_letter;
  composition_ = composition;
  mono_mass_ = mono_mass;
  avg_mass_ = avg_mass;
}

Acid::Acid (xercesc::DOMElement *element) {

  name_ = getChildValue(element, "name");
  one_letter_ = getChildValue(element, "one_letter");
  three_letter_ = getChildValue(element, "three_letter");
  composition_ = getChildValue(element, "composition");
  mono_mass_ = getDoubleChildValue(element, "mono_mass");
  avg_mass_ = getDoubleChildValue(element, "avg_mass");
}


}

int main(int argc, char** argv) {
  std::string value;
  proteomics::XmlDOMParser* parser = proteomics::getXmlDOMInstance();
  if (parser) {
    proteomics::XmlDOMDocument* doc 
        = new proteomics::XmlDOMDocument(parser, "./acid.xml");
    if (doc) {
      int acid_num = doc->getChildCount("amino_acid_list", 0, "amino_acid");
      for (int i = 0; i < acid_num; i++) {
        xercesc::DOMElement* element = doc->getElement("amino_acid", i);
        proteomics::Acid* acid = new proteomics::Acid(element);
        std::cout << acid->getName() << "\n";
        delete acid;
      }
      delete doc;
    }
    delete parser;
  }
  exit(0);
}

