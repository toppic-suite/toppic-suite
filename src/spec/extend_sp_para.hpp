#ifndef PROT_EXTEND_SP_PARA_HPP_
#define PROT_EXTEND_SP_PARA_HPP_

#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class ExtendSpPara {
 public:
	ExtendSpPara(double extend_min_mass, std::vector<double> ext_offsets);

	ExtendSpPara(xercesc::DOMElement* element);

  xercesc::DOMElement* toXml(XmlDOMDocument* xml_doc);

	// if the mass is smaller than extendMinMass, the peak is not extended 
	double extend_min_mass_;
  std::vector<double> ext_offsets_;
	
};

} /* name_space */

#endif
