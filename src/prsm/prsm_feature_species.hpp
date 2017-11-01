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


#ifndef PROT_PRSM_PRSM_FEATURE_SPECIES_HPP_
#define PROT_PRSM_PRSM_FEATURE_SPECIES_HPP_

#include <map>
#include <string>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

class PrsmFeatureSpecies {
 public:
  PrsmFeatureSpecies(const std::string &db_file_name,
                     const std::string &spec_file_name,
                     const std::string &feature_file_name,
                     const std::string &input_file_ext,
                     const std::string &output_file_ext,
                     const ModPtrVec &fix_mod_ptr_vec,
                     double prec_error_tole,
                     PrsmParaPtr prsm_para_ptr):
      db_file_name_(db_file_name),
      spec_file_name_(spec_file_name),
      feature_file_name_(feature_file_name),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext),
      fix_mod_ptr_vec_(fix_mod_ptr_vec),
      prec_error_tole_(prec_error_tole),
      prsm_para_ptr_(prsm_para_ptr) {}

  void process();

 private:
  std::string db_file_name_;
  std::string spec_file_name_;
  std::string feature_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  ModPtrVec fix_mod_ptr_vec_;
  double prec_error_tole_;
  PrsmParaPtr prsm_para_ptr_;

  void setProtId(PrsmStrPtrVec & prsm_ptrs);

  void setSpeciesId(PrsmStrPtrVec & prsm_ptrs);
};

typedef std::shared_ptr<PrsmFeatureSpecies> PrsmFeatureSpeciesPtr;

}  // namespace prot

#endif
