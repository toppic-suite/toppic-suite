#ifndef PROT_BASE_CHANGE_TYPE_HPP_
#define PROT_BASE_CHANGE_TYPE_HPP_

#include <memory>
#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class ChangeType;
typedef std::shared_ptr<ChangeType> ChangeTypePtr;

class ChangeType {
 public:
  static const ChangeTypePtr INPUT;
  static const ChangeTypePtr FIXED;
  static const ChangeTypePtr PROTEIN_VARIABLE;
  static const ChangeTypePtr VARIABLE;
  static const ChangeTypePtr UNEXPECTED;

  std::string getName() {return name_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static ChangeTypePtr getChangeTypePtrFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "change_type";}

 private:
  std::string name_;
  ChangeType(std::string name): name_(name) {};
};

typedef std::vector<ChangeTypePtr> ChangeTypePtrVec;

}

#endif

