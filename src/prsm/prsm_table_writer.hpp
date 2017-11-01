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


#ifndef PROT_PRSM_PRSM_TABLE_WRITER_HPP_
#define PROT_PRSM_PRSM_TABLE_WRITER_HPP_

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <map>

#include "base/string_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class PrsmTableWriter {
 public:
  PrsmTableWriter(PrsmParaPtr prsm_para_ptr, 
                  std::map<std::string, std::string> arguments,
                  const std::string &input_file_ext, 
                  const std::string &output_file_ext);
  void write();

  void writePrsm(std::ofstream &file, PrsmPtr prsm_ptr);

 private:
  PrsmParaPtr prsm_para_ptr_;
  std::map<std::string, std::string> arguments_;
  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<PrsmTableWriter> PrsmTableWriterPtr;

} /* namespace prot */

#endif /* TABLE_WRITER_HPP_ */
