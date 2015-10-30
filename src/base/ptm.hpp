/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROT_PTM_HPP_
#define PROT_PTM_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/acid.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Ptm;

class Ptm {
 public:
  Ptm(const std::string &name, const std::string &abbr_name,
      double mono_mass, int unimod_id, 
      const std::string &n_term_acid_str,
      const std::string &c_term_acid_str, 
      const std::string &anywhere_acid_str);

  const std::string& getName() {return name_;}

  const std::string& getAbbrName() {return abbr_name_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  const AcidPtrVec& getNTermAcids() {return n_term_acids_;}

  const AcidPtrVec& getCTermAcids() {return c_term_acids_;}

  const AcidPtrVec& getAnywhereAcids() {return anywhere_acids_;}

  int getUnimodId() {return unimod_id_;}

  static std::string getXmlElementName() {return "modification";}

  void appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  /* Full name */
  std::string name_;
  // abbrevation name 
  std::string abbr_name_;
  /* monoisotopic mass */
  double mono_mass_;
  // unimod id
  int unimod_id_;
  // possible acids
  AcidPtrVec n_term_acids_;
  AcidPtrVec c_term_acids_;
  AcidPtrVec anywhere_acids_;
};

typedef std::shared_ptr<Ptm> PtmPtr;
typedef std::vector<PtmPtr> PtmPtrVec;
typedef std::pair<PtmPtr, PtmPtr> PtmPair;
typedef std::vector<PtmPair> PtmPairVec;

}

#endif
