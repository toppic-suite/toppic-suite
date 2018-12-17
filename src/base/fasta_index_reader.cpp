//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#include <string>
#include <vector>

#include "base/logger.hpp"
#include "base/fasta_index_reader.hpp"

namespace toppic {

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
  if (seq_len < 0) {
    LOG_ERROR("Failed to fetch sequence " << name);
  }
  std::string ori_seq(seq);

  free(seq);

  return std::make_shared<FastaSeq>(name, desc, ori_seq);
}

std::vector<FastaSubSeqPtr> FastaIndexReader::readFastaSubSeqVec(const std::string & name,
                                                                 const std::string & desc) {
  int seq_len;
  char * seq = fai_fetch(fai_, name.c_str(), &seq_len);
  if (seq_len < 0) {
    LOG_ERROR("Failed to fetch sequence " << name);
  }

  std::string ori_seq(seq);
  free(seq);

  int N = 2000;
  std::vector<FastaSubSeqPtr> fasta_seq_vec;
  if (seq_len < N) {
    fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(name, desc, ori_seq, 0));
  } else {
    int k = seq_len / N;
    for (int i = 0; i <= k; i++) {
      std::string sub_seq;
      if (N * (i + 1) > seq_len) {
        sub_seq = ori_seq.substr(seq_len - N, N);
        fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(name, desc, sub_seq, seq_len - N));
      } else {
        sub_seq = ori_seq.substr(i * N, N);
        fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(name, desc, sub_seq, i * N));
      }
      if (i != k) {
        if (N * (i + 1.5) < seq_len) {
          sub_seq = ori_seq.substr((i + 0.5) * N, N);
          fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(name, desc, sub_seq, (i + 0.5) * N));
        }
      }
    }
  }
  return fasta_seq_vec;
}

}  // namespace toppic
