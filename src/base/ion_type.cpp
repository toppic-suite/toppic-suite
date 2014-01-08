#include "base/mass_constant.hpp"
#include "base/ion_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

IonType::IonType(std::string name, bool n_term, double shift) {
  name_ = name;
  n_term_ = n_term;
  shift_ = shift;
  if (n_term) {
	  b_y_shift_ = shift_;
  }
  else {
	  b_y_shift_ = shift_ - MassConstant::getWaterMass();
  }
}

void IonType::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
	xercesc::DOMElement* element = xml_doc->createElement("ion_type");
	xml_doc->addElement(element, "name", name_.c_str());
	std::string str = convertToString(n_term_);
	xml_doc->addElement(element, "n_term", str.c_str());
	str = convertToString(shift_);
	xml_doc->addElement(element, "shift", str.c_str());
	str = convertToString(b_y_shift_);
	xml_doc->addElement(element, "b_y_shift_", str.c_str());
	parent->appendChild(element);
}

IonTypePtrVec getIonTypePtrVecInstance(const std::string file_name){
	IonTypePtrVec ion_type_ptr_vec;
	prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int ion_type_num = getChildCount(parent, "ion_type");
    for (int i = 0; i < ion_type_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "ion_type", i);
      std::string name = getChildValue(element, "name", 0);
      bool n_term = getBoolChildValue(element, "n_term", 0);
      double shift = getDoubleChildValue(element, "shift", 0);
      ion_type_ptr_vec.push_back(IonTypePtr(new IonType(name, n_term, shift)));
    }
  }
	return ion_type_ptr_vec;
}

IonTypePtr getIonTypePtrByName(IonTypePtrVec &ion_type_list, const std::string &name){
	for (unsigned int i = 0; i < ion_type_list.size(); i++) {
	    std::string n = ion_type_list[i]->getName();
	    if (n == name) {
	      return ion_type_list[i];
	    }
	  }
	  return IonTypePtr(nullptr);
}

}
