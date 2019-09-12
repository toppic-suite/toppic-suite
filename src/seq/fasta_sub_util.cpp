//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "seq/fasta_sub_util.hpp" 

namespace toppic {

namespace fasta_sub_util {

std::vector<FastaSubSeqPtr> breakSeq(FastaSeqPtr seq_ptr) {
  int N = 2000;
  std::vector<FastaSubSeqPtr> fasta_seq_vec;
  int seq_len = seq_ptr->getAcidPtmPairLen();
  if (seq_len < N) {
    fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(seq_ptr, 0, seq_len));
  } 
  else {
    int k = seq_len / N;
    for (int i = 0; i <= k; i++) {
      if (N * (i + 1) > seq_len) {
        int sub_seq_len = seq_len - N * i;
        fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(seq_ptr, i*N, sub_seq_len));
      } else {
        fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(seq_ptr, i*N, N));
      }
      if (i != k) {
        if (N * (i + 1.5) < seq_len) {
          fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(seq_ptr, i*N + N/2, N));
        }
      }
    }
  }
  return fasta_seq_vec;
}

}  // namespace fasta_sub_util

}  // namespace toppic
