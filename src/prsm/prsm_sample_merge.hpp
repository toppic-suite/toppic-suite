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


#ifndef PROT_PRSM_PRSM_SAMPLE_MERGE_HPP_
#define PROT_PRSM_PRSM_SAMPLE_MERGE_HPP_

#include <map>
#include <string>
#include <vector>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

class PrsmSampleMerge {
 public:
  PrsmSampleMerge(const std::string &database_file_name,
                  const std::vector<std::string> &input_file_names,
                  const std::string &output_file_name,
                  const std::string &fix_mod,
                  double error_tole):
      database_file_name_(database_file_name),
      input_file_names_(input_file_names),
      output_file_name_(output_file_name),
      fix_mod_(fix_mod),
      error_tole_(error_tole) {}

  void process();

 private:
  std::string database_file_name_;
  std::vector<std::string> input_file_names_;
  std::string output_file_name_;
  std::string fix_mod_;
  double error_tole_;
};

typedef std::shared_ptr<PrsmSampleMerge> PrsmSampleMergePtr;

}  // namespace prot

#endif /* PRSM_SAMPLE_MERGE_HPP_ */
