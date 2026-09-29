// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_PARA_PRSM_PARA_HPP_
#define TOPPIC_PARA_PRSM_PARA_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "common/base/mod.hpp"
#include "common/base/prot_mod.hpp"
#include "common/base/proteoform_type.hpp"
#include "para/sp_para.hpp"

namespace toppic {

class PrsmPara {
 public:
  explicit PrsmPara(const std::map<std::string, std::string>& arguments);

  const std::string& getOriDbName() const { return ori_db_name_; }

  const std::string& getSearchDbFileNameWithFolder() const {
    return search_db_file_name_with_folder_;
  }

  std::string getDbIndexDir() const { return ori_db_name_ + "_idx"; }

  const std::string& getSpectrumFileName() const { return spec_file_name_; }

  const std::string& getResourceDir() const { return resource_dir_; }

  int getGroupSpecNum() const { return group_spec_num_; }

  const ModPtrVec& getFixModPtrVec() const { return fix_mod_list_; }

  const ProtModPtrVec& getProtModPtrVec() const { return prot_mod_list_; }

  const SpParaPtr& getSpParaPtr() const { return sp_para_ptr_; }

  bool allowProtType(const ProteoformTypePtr& type_ptr) const;

  bool allowCompleteProt() const {
    return allowProtType(ProteoformType::COMPLETE);
  }

  bool allowPrefixProt() const { return allowProtType(ProteoformType::PREFIX); }

  bool allowSuffixProt() const { return allowProtType(ProteoformType::SUFFIX); }

  bool allowInternalProt() const {
    return allowProtType(ProteoformType::INTERNAL);
  }

  // The multiple-unexpected-shift ("slow filter") search computes
  // COMPLETE/PREFIX/SUFFIX/INTERNAL for a candidate with one fused dynamic
  // program (search/ptmsearch/ps_align_multi_type.hpp) instead of one PsAlign
  // pass per type. Default true; toppic hidden --disable-multi-type-align
  // falls back to the original one-pass-per-type PsAlign.
  bool useMultiTypeAlign() const { return use_multi_type_align_; }

  // The number of unexpected mass shifts of TopMG's proteoform-graph
  // alignment, separate from toppic's --num-shift, which governs toppic's own
  // one-shift and multiple-shift stages. topmg sets it with its --num-shift
  // (default 1); toppic sets it with the hidden --topmg-num-shift of its
  // --topmg stage.
  int getTopmgShiftNum() const { return topmg_shift_num_; }

  void setTopmgShiftNum(int shift_num) { topmg_shift_num_ = shift_num; }

 private:
  std::string ori_db_name_;

  std::string search_db_file_name_;

  std::string search_db_file_name_with_folder_;

  std::string spec_file_name_;

  std::string resource_dir_;

  ModPtrVec fix_mod_list_;

  ProtModPtrVec prot_mod_list_;

  ProteoformTypePtrVec prot_type_list_;

  int group_spec_num_;

  /** spectrum parameters */
  SpParaPtr sp_para_ptr_;

  bool use_multi_type_align_ = true;

  int topmg_shift_num_ = 1;
};

using PrsmParaPtr = std::shared_ptr<PrsmPara>;
using PrsmParaPtrVec = std::vector<PrsmParaPtr>;

} /* namespace toppic */

#endif /* TOPPIC_PARA_PRSM_PARA_HPP_ */
