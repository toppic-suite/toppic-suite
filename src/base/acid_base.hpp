// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_BASE_ACID_BASE_HPP_
#define PROT_BASE_ACID_BASE_HPP_

#include "base/acid.hpp"

namespace prot {

class AcidBase {
 public:
  static void initBase(const std::string &file_name);

  static const AcidPtrVec& getBaseAcidPtrVec() {return acid_ptr_vec_;}

  static AcidPtr getEmptyAcidPtr() {return empty_acid_ptr_;}

  // Returns an amino acid based on the the name. Returns null if the amino
  // acid name does not exist.
  static AcidPtr getAcidPtrByName(const std::string &name);

  // Returns an amino acid based on the one letter representation. Returns
  // null if the one letter representation does not exist.
  static AcidPtr getAcidPtrByOneLetter(const std::string &one_letter);

  // Returns an amino acid based on the three letter representation. Returns
  // null if the three letter representation does not exist.
  static AcidPtr getAcidPtrByThreeLetter(const std::string &three_letter);

  // Checks if the list contains an amino acid with the specific name.
  static bool containsName(const std::string &name);

  // Checks if the list contains an amino acid with the specific one letter
  // representation.
  static bool containsOneLetter(const std::string &one_letter);

  // Checks if the list contains an amino acid with the specific three letter
  // representation.
  static bool containsThreeLetter(const std::string &three_letter);

  static AcidPtr getAcidPtrFromXml(xercesc::DOMElement * element);

 private:
  static AcidPtrVec acid_ptr_vec_;
  static AcidPtr empty_acid_ptr_;
};

}
#endif
