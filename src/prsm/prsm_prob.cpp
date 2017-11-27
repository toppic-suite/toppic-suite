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


#include <string>

#include "base/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_prob.hpp"

namespace prot {

void PrsmProb::process() {
  std::string input_file_name = FileUtil::basename(spec_file_name_)+ "." + input_file_ext_;
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name_);
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec_);

  PrsmXmlWriter all_writer(FileUtil::basename(spec_file_name_) + "." + output_file_ext_);

  while (prsm_ptr != nullptr) {
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
