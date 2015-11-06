#ifndef PROT_BASE_FASTA_INDEX_READER_HPP_
#define PROT_BASE_FASTA_INDEX_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "htslib/faidx.h"
#include "base/fasta_seq.hpp"
#include "base/string_util.hpp"

namespace prot {

class FastaIndexReader {
 public:
  FastaIndexReader(const std::string &file_name);
  ~FastaIndexReader();

  FastaSeqPtr readFastaSeq(const std::string &name,
                           const std::string &desc);

 private:
  faidx_t *fai_;
};

typedef std::shared_ptr<FastaIndexReader> FastaIndexReaderPtr;

}  //namepace prot

#endif
