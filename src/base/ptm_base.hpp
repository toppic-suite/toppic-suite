//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


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

  static bool isEmptyPtmPtr(PtmPtr ptm_ptr) {return ptm_ptr == empty_ptm_ptr_;}

  static PtmPtr getPtmPtr_Acetylation() {return acetylation_ptr_;}
  static PtmPtr getPtmPtr_C57() {return c57_ptr_;}
  static PtmPtr getPtmPtr_C58() {return c58_ptr_;}
  /**
   * Returns a PTM based on the abbreviation name. Returns null if the
   * abbreviation name does not exist.
   */
  static PtmPtr getPtmPtrByAbbrName(const std::string &abbr_name);

  static PtmPtr getPtmPtr(PtmPtr p);

  /**
   * Checks if the list contains an amino acid with the specific name.
   */
  static bool containsAbbrName(const std::string &abbr_name);

  static PtmPtr getPtmPtrFromXml(xercesc::DOMElement * element);

 private:
  static PtmPtrVec ptm_ptr_vec_;
  static PtmPtr empty_ptm_ptr_;
  static PtmPtr acetylation_ptr_;
  static PtmPtr c57_ptr_;
  static PtmPtr c58_ptr_;

  static std::string getAcetylationAbbrName() {return "Acetyl";}
  static std::string getC57AbbrName() {return "Carbamidomethylation";}
  static std::string getC58AbbrName() {return "Carboxymethyl";}
};

}

#endif

