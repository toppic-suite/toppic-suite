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


#ifndef PROT_PROTEO_ANNO_HPP_
#define PROT_PROTEO_ANNO_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>

#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

#include "prsm/prsm_para.hpp"
#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"

namespace prot {

class ProteoAnno {

 public:
  ProteoAnno(const ModPtrVec &fix_mod_ptr_vec, 
             const ProtModPtrVec &prot_mod_ptr_vec, 
             const ModPtrVec &var_mod_ptr_vec);

  void anno(const std::string &seq, bool is_complete = true);

  const ResiduePtrVec2D& getResiduePtrVec() {return res_vec_2d_;}
  const std::vector<std::vector<int>>& getChangeVec() {return change_vec_2d_;}

  const ResiduePtrVec& getResiduePtrVec(int i) {return res_vec_2d_[i];}
  const std::vector<int>& getChangeVec(int i) {return change_vec_2d_[i];}

  int getLen() {return res_vec_2d_.size();}
  bool isNme() {return is_nme_;}

 private:
  ResiduePtrVec2D res_vec_2d_;
  std::vector<std::vector<int>> change_vec_2d_;

  ModPtrVec fix_mod_ptr_vec_;
  ProtModPtrVec prot_mod_ptr_vec_;
  ModPtrVec var_mod_ptr_vec_;

  std::map<AcidPtr, ResiduePtrVec> ptm_map_;
  bool is_nme_;

};

typedef std::shared_ptr<ProteoAnno> ProteoAnnoPtr;

}

#endif
