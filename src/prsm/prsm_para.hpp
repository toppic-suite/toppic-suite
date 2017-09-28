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


#ifndef PROT_PRSM_PRSM_PARA_HPP_
#define PROT_PRSM_PRSM_PARA_HPP_

#include <string>
#include <map>

#include "base/mod.hpp"
#include "base/prot_mod.hpp"
#include "base/activation.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class PrsmPara {
 public:
  PrsmPara(std::map<std::string,std::string> &arguments);
  std::string getSearchDbFileName() {return search_db_file_name_;}
  std::string getSpectrumFileName() {return spec_file_name_;}
  std::string getExeDir(){return exe_dir_;}
  int getErrorTolerance(){return errorTolerance_;}
  int getGroupSpecNum() {return group_spec_num_;}
  const ModPtrVec& getFixModPtrVec() {return fix_mod_list_;}
  const ProtModPtrVec& getProtModPtrVec() {return prot_mod_list_;}
  SpParaPtr getSpParaPtr() {return sp_para_ptr_;}
  bool doLocaliztion() {return localization_;}

 private:
  std::string search_db_file_name_;
  std::string spec_file_name_;
  std::string exe_dir_;
  int errorTolerance_;

  ModPtrVec fix_mod_list_;
  ProtModPtrVec prot_mod_list_;

  int group_spec_num_;

  bool localization_;

  /** spectrum parameters */
  SpParaPtr sp_para_ptr_;
};

typedef std::shared_ptr<PrsmPara> PrsmParaPtr;

} /* namespace prot */

#endif /* PRSM_PARA_HPP_ */
