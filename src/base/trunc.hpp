#ifndef PROT_BASE_TRUNC_HPP_
#define PROT_BASE_TRUNC_HPP_

#include <string>
#include "base/residue.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Trunc {
 public:
  Trunc(const std::string &name, int trunc_len, 
        const std::string &trunc_residues,
        const std::string &allow_first_remain_residues_);

  Trunc(xercesc::DOMElement* element); 

  const std::string& getName() {return name_;}

  int getTruncLen() {return trunc_len_;}

  const ResiduePtrVec& getTruncResiduePtrVec() {return trunc_residue_ptr_vec_;}

  const ResiduePtrVec& getAllowFirstRemainResiduePtrs() {return allow_first_remain_residue_ptrs_;}

  double getShift() {return shift_;}

  static std::string getNameFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "truncation";}

 private:
  std::string name_;
  int trunc_len_;
  ResiduePtrVec trunc_residue_ptr_vec_;
  ResiduePtrVec allow_first_remain_residue_ptrs_;
  double shift_;
};

typedef std::shared_ptr<Trunc> TruncPtr;
typedef std::vector<TruncPtr> TruncPtrVec;

}

#endif
