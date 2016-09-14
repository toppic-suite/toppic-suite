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


#ifndef PROT_BASE_MOD_BASE_HPP_
#define PROT_BASE_MOD_BASE_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/mod.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class ModBase {
 public:
  static void initBase(const std::string &file_name);

  static const ModPtrVec& getBaseModPtrVec() {return mod_ptr_vec_;}

  static ModPtr getNoneModPtr() {return none_mod_ptr_;}

  static ModPtr getC57ModPtr() {return c57_mod_ptr_;}

  static ModPtr getC58ModPtr() {return c58_mod_ptr_;}

  static ModPtr getBaseModPtr(ModPtr mod_ptr);

  static ModPtr getBaseModPtr(ResiduePtr ori_residue, ResiduePtr mod_residue);

  static bool isNoneModPtr(ModPtr mod_ptr) {return mod_ptr == none_mod_ptr_;}

  static ModPtr getModPtrFromXml(xercesc::DOMElement * element);

 private:
  static ModPtrVec mod_ptr_vec_;
  static ModPtr none_mod_ptr_;
  //C57
  static ModPtr c57_mod_ptr_;
  //C58
  static ModPtr c58_mod_ptr_;
};

}

#endif

