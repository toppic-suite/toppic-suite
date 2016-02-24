#ifndef PROT_BASE_FASTA_READER_HPP_
#define PROT_BASE_FASTA_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "base/fasta_seq.hpp"
#include "base/string_util.hpp"

namespace prot {

class FastaReader {
 public:
  FastaReader(const std::string &file_name);

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
