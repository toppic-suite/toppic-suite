#ifndef PROT_EXTREME_VALUE_HPP_
#define PROT_EXTREME_VALUE_HPP_

#include <memory>
#include <vector>
#include "base/xml_dom_document.hpp"

namespace prot {

class ExtremeValue {
 public:
  ExtremeValue (double one_match_prob, double test_num, double adjust_factor);

  ExtremeValue (xercesc::DOMElement* element);

    double getPValue() {return p_value_;}
    
    double getEValue() {return e_value_;}

    double getOneProtProb() {return one_match_prob_;}

    double getTestNum() {return test_num_;}

    double getAdjustFactor() { return adjust_factor_;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
    /* oneMatchProb considers expected PTMs, but not unexpected PTMs */
    double one_match_prob_;
    double test_num_;
    double adjust_factor_;
    double p_value_;
    double e_value_;

  void init();
};

typedef std::shared_ptr<ExtremeValue> ExtremeValuePtr;
typedef std::vector<ExtremeValuePtr> ExtremeValuePtrVec;

}

#endif
