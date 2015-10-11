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
  ProteoAnno(const ResiduePtrVec &fix_mod_res_ptr_vec, 
             const ProtModPtrVec &prot_mod_ptr_vec, 
             const ResiduePtrVec &residue_mod_ptr_vec);

  void anno(const std::string &seq);

  const ResiduePtrVec2D& getResiduePtrVec() {return res_vec_2d_;}
  const std::vector<std::vector<int>>& getChangeVec() {return change_vec_2d_;}

  const ResiduePtrVec& getResiduePtrVec(int i) {return res_vec_2d_[i];}
  const std::vector<int>& getChangeVec(int i) {return change_vec_2d_[i];}

  int getLen() {return res_vec_2d_.size();}
  bool isNme() {return is_nme_;}

 private:
  ResiduePtrVec2D res_vec_2d_;
  std::vector<std::vector<int>> change_vec_2d_;

  ResiduePtrVec fix_mod_res_ptr_vec_;
  ProtModPtrVec prot_mod_ptr_vec_;
  ResiduePtrVec residue_mod_ptr_vec_;

  std::map<AcidPtr, ResiduePtrVec> ptm_map_;
  bool is_nme_;

};

typedef std::shared_ptr<ProteoAnno> ProteoAnnoPtr;

}

#endif
