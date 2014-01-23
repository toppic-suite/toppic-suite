#ifndef PROT_ION_TYPE_HPP_
#define PROT_ION_TYPE_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/xml_dom_document.hpp"

namespace prot {

class IonType;
typedef std::shared_ptr<IonType> IonTypePtr;

class IonType {
 public: 
  IonType(std::string name, bool n_term, double shift);

  std::string getName() {return name_;}

  bool isNTerm() {return n_term_;}

  double getShift() {return shift_;}

  double getBYShift() {return b_y_shift_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);
 private:
  // ion name
  std::string name_;
  // A B C are n-terminal ions and X Y Z are non-n-terminal ions 
  bool n_term_;
  /**
   * Shift stands for the shift of the ion compared to M. For example, the
   * shift for B ion is 0, and the shift for Y ion is 18 (chrg 0);
   */
  double shift_;

  // shifts compared to b or y-ions
  double b_y_shift_;
};

typedef std::shared_ptr<IonType> IonTypePtr;
typedef std::vector<IonTypePtr> IonTypePtrVec;

IonTypePtrVec getIonTypePtrVecInstance(const std::string file_name);
IonTypePtr getIonTypePtrByName(IonTypePtrVec &ion_type_list, const std::string &name);

}

#endif
