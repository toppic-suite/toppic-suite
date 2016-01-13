/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROT_BASE_PTM_HPP_
#define PROT_BASE_PTM_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/acid.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Ptm;
typedef std::shared_ptr<Ptm> PtmPtr;

class Ptm {
 public:
  Ptm(const std::string &name, const std::string &abbr_name,
      double mono_mass, int unimod_id, 
      const std::string &n_term_residue_str,
      const std::string &c_term_residue_str, 
      const std::string &anywhere_residue_str);

  Ptm(xercesc::DOMElement* element); 

  const std::string& getName() {return name_;}

  const std::string& getAbbrName() {return abbr_name_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  int getUnimodId() {return unimod_id_;}

  void appendAbbrNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getAbbrNameFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "ptm";}

  // comparison function
  static bool cmpMassInc(const PtmPtr &a, const PtmPtr &b) {
    return a->getMonoMass() < b->getMonoMass();
  }

  bool isSame(PtmPtr ptm_ptr) { return abbr_name_ == ptm_ptr->getAbbrName();}

 private:
  /* Full name */
  std::string name_;
  // abbrevation name 
  std::string abbr_name_;
  /* monoisotopic mass */
  double mono_mass_;
  // unimod id
  int unimod_id_;
};

typedef std::vector<PtmPtr> PtmPtrVec;

}

#endif
