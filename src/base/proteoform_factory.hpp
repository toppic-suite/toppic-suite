// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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
