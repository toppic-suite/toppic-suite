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


#include "base/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_prob.hpp"

namespace prot {

PrsmProb::PrsmProb(const std::string &db_file_name, 
                   const std::string &spec_file_name, 
                   const ModPtrVec &fix_mod_ptr_vec,
                   const std::string &in_file_ext,
                   const std::string &out_file_ext, 
                   double K1, double K2, double pref, double inte): 
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    fix_mod_ptr_vec_(fix_mod_ptr_vec),
    input_file_ext_(in_file_ext),
    output_file_ext_(out_file_ext),
    K1_(K1),
    K2_ (K2),
    pref_(pref),
    inte_(inte) {
    }

void PrsmProb::process() {

  std::string input_file_name = FileUtil::basename(spec_file_name_)+ "." + input_file_ext_;
  FastaIndexReaderPtr seq_reader(new FastaIndexReader(db_file_name_));
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec_);

  PrsmXmlWriter all_writer(FileUtil::basename(spec_file_name_)+"."+output_file_ext_);

  while(prsm_ptr != nullptr) {
    int shift_num = prsm_ptr->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED);
    AlignTypePtr type_ptr = prsm_ptr->getProteoformPtr()->getAlignType();
    ExtremeValuePtr prob_ptr = prsm_ptr->getExtremeValuePtr();
    if (shift_num == 1) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * K1_);
    }
    if (shift_num == 2) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * K2_);
    }
    if (type_ptr == AlignType::PREFIX || type_ptr == AlignType::SUFFIX) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * pref_);
    }

    if (type_ptr == AlignType::INTERNAL) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * inte_);
    }
    all_writer.write(prsm_ptr);
    prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec_);
  }

  prsm_reader.close();
  all_writer.close();
}

} /* namespace prot */
