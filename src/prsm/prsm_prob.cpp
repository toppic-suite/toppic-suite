//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "common/util/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_prob.hpp"

namespace toppic {

PrsmProb::PrsmProb(const std::string &db_file_name,
                   const std::string &spec_file_name,
                   const ModPtrVec &fix_mod_ptr_vec,
                   const std::string &in_file_ext,
                   const std::string &out_file_ext,
                   double K1, double K2,
                   double pref, double inte):
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    fix_mod_ptr_vec_(fix_mod_ptr_vec),
    input_file_ext_(in_file_ext),
    output_file_ext_(out_file_ext),
    K1_(K1),
    K2_(K2),
    pref_(pref),
    inte_(inte) {}

void PrsmProb::process() {
  std::string input_file_name = file_util::basename(spec_file_name_)+ "." + input_file_ext_;
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name_);
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec_);

  PrsmXmlWriter all_writer(file_util::basename(spec_file_name_) + "." + output_file_ext_);

  while (prsm_ptr != nullptr) {
    int shift_num = prsm_ptr->getProteoformPtr()->getMassShiftNum(AlterType::UNEXPECTED);
    ProteoformTypePtr type_ptr = prsm_ptr->getProteoformPtr()->getProteoformType();
    ExtremeValuePtr prob_ptr = prsm_ptr->getExtremeValuePtr();
    if (shift_num == 1) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * K1_);
    }
    if (shift_num == 2) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * K2_);
    }
    if (type_ptr == ProteoformType::PREFIX || type_ptr == ProteoformType::SUFFIX) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * pref_);
    }

    if (type_ptr == ProteoformType::INTERNAL) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * inte_);
    }
    all_writer.write(prsm_ptr);
    prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec_);
  }

  prsm_reader.close();
  all_writer.close();
}

} /* namespace toppic */
