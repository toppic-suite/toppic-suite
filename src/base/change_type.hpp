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

  int getId() {return id_;}
  std::string getName() {return name_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static ChangeTypePtr getChangeTypePtrFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "change_type";}

 private:
  int id_;
  std::string name_;
  ChangeType(int id, std::string name): id_(id), name_(name) {};
};

typedef std::vector<ChangeTypePtr> ChangeTypePtrVec;

}

#endif

