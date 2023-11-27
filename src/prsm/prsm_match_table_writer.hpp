//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_PRSM_PRSM_MATCH_TABLE_WRITER_HPP_
#define TOPPIC_PRSM_PRSM_MATCH_TABLE_WRITER_HPP_

#include "para/prsm_para.hpp"
#include "prsm/prsm.hpp"
#include "prsm/search_fasta_match.hpp"

namespace toppic {

class PrsmMatchTableWriter {
 public:
  PrsmMatchTableWriter(PrsmParaPtr prsm_para_ptr,  
                       std::string argu_str,
                       const std::string &input_file_ext, 
                       const std::string &output_file_ext,
                       bool write_multiple_matches);

  void write();

  void writePrsm(std::ofstream &file, PrsmPtr prsm_ptr);
  void writePrsmStandardFormat(std::ofstream &file, PrsmPtr prsm_ptr);
  std::string formatSeq(std::string seq);

  void setOutputName(std::string output_file_ext) {output_file_ext_ = output_file_ext;}

  void setWriteMultiMatches(bool write_multiple_matches) {write_multiple_matches_ = write_multiple_matches;}

  void setOutputName(std::string output_file_ext) {output_file_ext_ = output_file_ext;}

  void setWriteMultiMatches(bool write_multiple_matches) {write_multiple_matches_ = write_multiple_matches;}

 private:
  PrsmParaPtr prsm_para_ptr_;

  std::string input_file_ext_;

  std::string argu_str_;

  std::string output_file_ext_;

  SearchFastaMatchPtr search_match_ptr_;

  bool write_multiple_matches_;
};

typedef std::shared_ptr<PrsmMatchTableWriter> PrsmMatchTableWriterPtr;

} /* namespace toppic */

#endif /* TABLE_WRITER_HPP_ */
