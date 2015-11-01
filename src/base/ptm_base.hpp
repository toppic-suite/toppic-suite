/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROT_BASE_PTM_BASE_HPP_
#define PROT_BASE_PTM_BASE_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class PtmBase {
 public:
  static void initBase(const std::string &file_name);

  static const PtmPtrVec& getBasePtmPtrVec() {return ptm_ptr_vec_;}

  static PtmPtr getEmptyPtmPtr() {return empty_ptm_ptr_;}

  static PtmPtr getPtmPtr_Acetylation() {return acetylation_ptr_;}
  /**
   * Returns a PTM based on the abbreviation name. Returns null if the
   * abbreviation name does not exist.
   */
  static PtmPtr getPtmPtrByAbbrName(const std::string &abbr_name);
  /**
   * Checks if the list contains an amino acid with the specific name.
   */
  static bool containsAbbrName(const std::string &abbr_name);

  //static PtmPtr addBasePtm(const std::string &abbr_name, double mono_mass);
  
  static PtmPtr getPtmPtrFromXml(xercesc::DOMElement * element);

 private:
  static PtmPtrVec ptm_ptr_vec_;
  static PtmPtr empty_ptm_ptr_;
  static PtmPtr acetylation_ptr_;

  static std::string getAcetylationAbbrName() {return "Acetylation";}
};

}

#endif

