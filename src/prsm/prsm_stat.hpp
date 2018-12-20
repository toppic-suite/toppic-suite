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


#ifndef PROT_PRSM_PRSM_STAT_HPP_
#define PROT_PRSM_PRSM_STAT_HPP_

#include <string>

#include "util/str_util.hpp"
#include "seq/proteoform.hpp"
#include "seq/fasta_reader.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm.hpp"

namespace toppic {

class PrsmStat {
 public:
  PrsmStat(PrsmParaPtr prsm_para_ptr, 
           const std::string &input_file_ext, 
           const std::string &output_file_ext);
  void process();

  void writePrsm(std::ofstream &file, PrsmPtr prsm_ptr);

 private:
  PrsmParaPtr prsm_para_ptr_;
  double min_mass_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  AminoAcidPtrVec acid_ptr_vec_;
};

typedef std::shared_ptr<PrsmStat> PrsmStatPtr;


} /* namespace toppic */

#endif /* PRSM_STAT_HPP_ */
