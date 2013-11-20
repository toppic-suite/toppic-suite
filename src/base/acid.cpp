/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */

#include "acid.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

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


