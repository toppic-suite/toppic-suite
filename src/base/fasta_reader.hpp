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
  /**
   * Constructs an instance with a File.
   **/
  FastaReader(const std::string &file_name);

  void setSeqId(int seq_id) {seq_id_ = seq_id;}

  /**
   * Read FASTA file and return next protein
   * name and sequence. 
   **/                  
  FastaSeqPtr getNextSeq();

 private:
  std::ifstream input_;
  std::string ori_name_;
  int seq_id_ = 0;
};

typedef std::shared_ptr<FastaReader> FastaReaderPtr;

}  //namepace prot

#endif
