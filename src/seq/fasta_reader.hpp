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

#ifndef TOPPIC_SEQ_FASTA_READER_HPP_
#define TOPPIC_SEQ_FASTA_READER_HPP_

#include <fstream>

#include "seq/fasta_seq.hpp"

namespace toppic {

class FastaReader {
 public:
  FastaReader(const std::string &file_name);

  ~FastaReader();

  // Read FASTA file and return next protein
  // name and sequence. 
  FastaSeqPtr getNextSeq();

  void close();

 private:
  std::ifstream input_;
  std::string ori_name_;
};

typedef std::shared_ptr<FastaReader> FastaReaderPtr;

}  //namepace prot

#endif
