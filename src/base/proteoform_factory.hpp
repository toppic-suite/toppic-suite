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


#ifndef PROT_BASE_PROTEOFORM_FACTORY_HPP_
#define PROT_BASE_PROTEOFORM_FACTORY_HPP_

#include "base/proteoform.hpp"
#include "base/fasta_index_reader.hpp"

namespace prot {

class ProteoformFactory {
 public:
  /* get db proteoform */
  static ProteoformPtr geneDbProteoformPtr(FastaSeqPtr seq_ptr, ModPtrVec fix_mod_list);

  /* generate a proteoform with protein mod */
  static ProteoformPtr geneProtModProteoform(ProteoformPtr db_form_ptr,
                                             ProtModPtr prot_mod_ptr);

  static ProteoformPtrVec geneProtModProteoform(ProteoformPtr db_form_ptr,
                                                const ProtModPtrVec &prot_mod_ptrs);

  static ProteoformPtrVec2D gene2DProtModProteoform(const ProteoformPtrVec &db_form_ptrs,
                                                    const ProtModPtrVec &prot_mod_ptrs);
  /*
   * get subproteoform. local_start and local_end are relatively to
   * the start position in the original proteoform
   */
  static ProteoformPtr geneSubProteoform(ProteoformPtr proteoform_ptr,
                                         int local_start, int local_end);

  /* generate a proteoform vector with protein mod */
  static ProteoformPtrVec geneProtModProteoform(const ProteoformPtrVec &ori_forms,
                                                const ProtModPtrVec &prot_mods);

  static ProteoformPtrVec readFastaToProteoformPtrVec(const std::string &file_name, 
                                                      const ModPtrVec &fix_mod_list);

  static ProteoformPtr readFastaToProteoformPtr(FastaIndexReaderPtr reader_ptr, 
                                                const std::string &seq_name,
                                                const std::string &seq_desc,
                                                const ModPtrVec &fix_mod_list);

};

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
