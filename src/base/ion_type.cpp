#include "mass_constant.hpp"
#include "ion_type.hpp"
#include "xml_dom_document.hpp"

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

IonType::IonType(xercesc::DOMElement *element) {

	name_ = getChildValue(element, "name");
	n_term_ = getBoolChildValue(element, "n_term");
	shift_ = getDoubleChildValue(element, "shift");
	if (n_term_) {
		b_y_shift_ = shift_;
	} else {
		b_y_shift_ = shift_ - MassConstant::getWaterMass();
	}
}

IonTypePtrVec getIonTypePtrVecInstance(const char* file_name){
	IonTypePtrVec ionType_ptr_vec;
	prot::XmlDOMParser* parser = prot::getXmlDOMInstance();
	if (parser) {
	    prot::XmlDOMDocument* doc = new prot::XmlDOMDocument(parser, file_name);
	    if (doc) {
	      int acid_num = doc->getChildCount("ion_type_list", 0, "ion_type");
	      for (int i = 0; i < acid_num; i++) {
	        xercesc::DOMElement* element = doc->getElement("ion_type", i);
	        ionType_ptr_vec.push_back(IonTypePtr(new IonType(element)));

	      }
	      delete doc;
	    }
	    delete parser;
	  }
	return ionType_ptr_vec;
}

IonTypePtr getIonTypePtrByName(IonTypePtrVec &ion_type_list, const std::string &name){
	for (unsigned int i = 0; i < ion_type_list.size(); i++) {
	    std::string n = ion_type_list[i]->getName();
	    if (n.compare(name) == 0) {
	      return ion_type_list[i];
	    }
	  }
	  return IonTypePtr(nullptr);
}

}
