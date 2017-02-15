// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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

