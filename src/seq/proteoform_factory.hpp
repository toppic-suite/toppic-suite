// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_SEQ_PROTEOFORM_FACTORY_HPP_
#define TOPPIC_SEQ_PROTEOFORM_FACTORY_HPP_

#include <string>

#include "seq/fasta_index_reader.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

namespace proteoform_factory {

ProteoformPtr geneDbProteoformPtr(const FastaSeqPtr& seq_ptr,
                                  const ModPtrVec& fix_mod_list);

ProteoformPtr geneDbProteoformPtr(const FastaSeqPtr& seq_ptr,
                                  const ModPtrVec& fix_mod_list, int start_pos);

// Generate a proteoform with protein N-terminal modification
ProteoformPtr geneProtModProteoform(const ProteoformPtr& db_form_ptr,
                                    const ProtModPtr& prot_mod_ptr);

ProteoformPtrVec geneProtModProteoform(const ProteoformPtr& db_form_ptr,
                                       const ProtModPtrVec& prot_mod_ptrs);

ProteoformPtrVec geneProtModProteoform(const ProteoformPtrVec& ori_forms,
                                       const ProtModPtrVec& prot_mods);

ProteoformPtrVec2D gene2DProtModProteoform(const ProteoformPtrVec& db_form_ptrs,
                                           const ProtModPtrVec& prot_mod_ptrs);

// Generate a proteoform with a truncation N-terminal modification: a copy of
// form_ptr whose N-terminus (the alpha-amine of its first residue) carries
// mod_ptr, a variable PTM of the proteoform N-terminus (Position N-term in a
// modification file). This is different from the protein N-terminal forms of
// geneProtModProteoform: the proteoform must have no protein N-terminal form
// and start at the third residue of the protein or later, and the first
// residue's side chain stays free for a variable PTM. Returns null when the
// modification cannot be added (Proteoform::setTruncNTermModPtr).
ProteoformPtr geneTruncNTermModProteoform(const ProteoformPtr& form_ptr,
                                          const ModPtr& mod_ptr);

// The proteoforms of geneTruncNTermModProteoform for the modifications of
// mod_ptrs that can be added
ProteoformPtrVec geneTruncNTermModProteoform(const ProteoformPtr& form_ptr,
                                             const ModPtrVec& mod_ptrs);

// Get subproteoform. local_start and local_end are relatively to
// the start position in the original proteoform. A subproteoform that
// starts at the first residue keeps the truncation N-terminal modification.
ProteoformPtr geneSubProteoform(const ProteoformPtr& proteoform_ptr,
                                const FastaSeqPtr& fasta_seq_ptr,
                                int local_start, int local_end);

ProteoformPtrVec readFastaToProteoformPtrVec(const std::string& file_name,
                                             const ModPtrVec& fix_mod_list);

ProteoformPtr readFastaToProteoformPtr(const FastaIndexReaderPtr& reader_ptr,
                                       const std::string& seq_name,
                                       const std::string& seq_desc,
                                       const ModPtrVec& fix_mod_list);

}  // namespace proteoform_factory

}  // namespace toppic

#endif
