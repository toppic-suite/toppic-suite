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


#ifndef PROT_BASE_FASTA_INDEX_READER_HPP_
#define PROT_BASE_FASTA_INDEX_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "htslib/faidx.h"
#include "base/fasta_seq.hpp"
#include "base/fasta_sub_seq.hpp"
#include "base/string_util.hpp"

namespace prot {

class FastaIndexReader {
 public:
  explicit FastaIndexReader(const std::string &file_name);
  
  ~FastaIndexReader();

  FastaSeqPtr readFastaSeq(const std::string &name,
                           const std::string &desc);

  std::vector<FastaSubSeqPtr> readFastaSubSeqVec(const std::string & name,
                                                 const std::string & desc);

 private:
  faidx_t *fai_;
};

typedef std::shared_ptr<FastaIndexReader> FastaIndexReaderPtr;

}  // namespace prot

#endif
