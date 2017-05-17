// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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
  return  std::make_shared<FastaSeq>(name, desc, ori_seq);
}

std::vector<FastaSeqPtr> FastaIndexReader::readFastaSeqVec(const std::string & name, 
                                                           const std::string & desc) {
  int seq_len;
  char * seq = fai_fetch(fai_, name.c_str(), &seq_len);
  if ( seq_len < 0 ) {
    LOG_ERROR("Failed to fetch sequence " << name);
  }

  std::string ori_seq(seq);
  free(seq);

  int N = 2000;
  std::vector<FastaSeqPtr> fasta_seq_vec;
  if (seq_len < N) {
    fasta_seq_vec.push_back(std::make_shared<FastaSeq>(name, desc, ori_seq));
  } else {
    int k = seq_len / N; 
    for (int i = 0; i <= k; i++) {
      std::string sub_seq;
      if (N * (i + 1) > seq_len) {
        sub_seq = ori_seq.substr(seq_len - N, N);
        fasta_seq_vec.push_back(std::make_shared<FastaSeq>(name, desc, sub_seq, seq_len - N));
      } else {
        sub_seq = ori_seq.substr(i * N, N);
        fasta_seq_vec.push_back(std::make_shared<FastaSeq>(name, desc, sub_seq, i * N));
      }
      if (i != k) {
        if (N * (i + 1.5) < seq_len) {
          sub_seq = ori_seq.substr((i + 0.5) * N, N);
          fasta_seq_vec.push_back(std::make_shared<FastaSeq>(name, desc, sub_seq, (i + 0.5) * N));
        }
      }
    }
  }
  return fasta_seq_vec;
}

}
