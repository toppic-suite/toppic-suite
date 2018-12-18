//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include "seq/fasta_reader.hpp"
#include "util/file_util.hpp"

namespace toppic {

class ProteoAnno {
 public:
  ProteoAnno(const ModPtrVec &fix_mod_ptr_vec,
             const ProtModPtrVec &prot_mod_ptr_vec,
             const ModPtrVec &var_mod_ptr_vec);

  void anno(const std::string &seq, bool is_complete = true);

  const ResiduePtrVec2D& getResiduePtrVec() {return res_vec_2d_;}

  const std::vector<std::vector<int> >& getChangeVec() {return shift_vec_2d_;}

  const ResiduePtrVec& getResiduePtrVec(int i) {return res_vec_2d_[i];}

  const std::vector<int>& getChangeVec(int i) {return shift_vec_2d_[i];}

  int getLen() {return res_vec_2d_.size();}

  bool isNme() {return is_nme_;}

 private:
  ResiduePtrVec2D res_vec_2d_;

  std::vector<std::vector<int> > shift_vec_2d_;

  ModPtrVec fix_mod_ptr_vec_;

  ProtModPtrVec prot_mod_ptr_vec_;

  ModPtrVec var_mod_ptr_vec_;

  std::map<AminoAcidPtr, ResiduePtrVec> ptm_map_;

  bool is_nme_;
};

typedef std::shared_ptr<ProteoAnno> ProteoAnnoPtr;

}  // namespace toppic

#endif
