#include "base/logger.hpp"
#include "base/fasta_index_reader.hpp"

namespace prot {

FastaIndexReader::FastaIndexReader(const std::string &file_name) {
  fai_ = fai_load(file_name.c_str());
}

FastaIndexReader::~FastaIndexReader() {
  fai_destroy(fai_);
}

FastaSeqPtr FastaIndexReader::readFastaSeq(const std::string &name, 
                                           const std::string &desc) {
  int seq_len;
  char *seq = fai_fetch(fai_, name.c_str(), &seq_len);
  if ( seq_len < 0 ) {
    LOG_ERROR("Failed to fetch sequence " << name);
  }
  std::string ori_seq(seq);
  free(seq);
  FastaSeqPtr fasta_seq_ptr(new FastaSeq(name, desc, ori_seq));
  return fasta_seq_ptr;
}

}
