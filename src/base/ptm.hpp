/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROT_PTM_HPP_
#define PROT_PTM_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Ptm;
typedef std::shared_ptr<Ptm> PtmPtr;
typedef std::vector<PtmPtr> PtmPtrVec;

class Ptm {
 public:
  Ptm(const std::string &abbr_name, 
      double mono_mass);

  /* Get amino acid composition. */
  std::string getAbbrName() { return abbr_name_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  /* Is it an empty PTM. */
  bool isEmpty();

  void appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  /* Abbreviation name of a PTM */
  std::string abbr_name_;
  /* monoisotopic mass */
  double mono_mass_;
};

/**
 * Returns a PTM based on the abbreviation name. Returns null if the
 * abbreviation name does not exist.
 */
PtmPtr getPtmPtrByAbbrName(PtmPtrVec &ptm_list, 
                           const std::string &abbr_name);

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containAbbrName(PtmPtrVec &ptm_list, const std::string &abbr_name);


PtmPtr addPtm(PtmPtrVec &ptm_list, std::string abbr_name,
              double mono_mass);

/* ptm factory */
class PtmFactory {
 public:
  static void initFactory(const std::string &file_name);
  static PtmPtrVec& getBasePtmPtrVec() {return ptm_ptr_vec_;}
  static PtmPtr findEmptyPtmPtr();
  static PtmPtr getBasePtmPtrByAbbrName(const std::string &abbr_name) {
    return getPtmPtrByAbbrName(ptm_ptr_vec_, abbr_name);
  }
  static bool baseContainAbbrName(const std::string &abbr_name) {
    return containAbbrName(ptm_ptr_vec_, abbr_name);
  }

 private:
  static PtmPtrVec ptm_ptr_vec_;
};

}
#endif

