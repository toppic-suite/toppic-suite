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

#include "seq/fasta_sub_util.hpp" 

namespace toppic {

namespace fasta_sub_util {

FastaSubSeqPtrVec breakSeq(FastaSeqPtr seq_ptr) {
  int N = 2000;
  return breakSeq(seq_ptr, N);
}

FastaSubSeqPtrVec breakSeq(FastaSeqPtr seq_ptr, int N) {
  std::vector<FastaSubSeqPtr> fasta_seq_vec;
  int seq_last_pos = seq_ptr->getAcidPtmPairLen() - 1;
  int sub_bgn_pos = 0;
  int sub_end_pos = N-1;
  bool last = false;
  do {
    int sub_len = N;
    if (seq_last_pos <= sub_end_pos) {
      sub_len = seq_last_pos - sub_bgn_pos + 1;
      last = true;
    }
    fasta_seq_vec.push_back(std::make_shared<FastaSubSeq>(seq_ptr, sub_bgn_pos, sub_len));
    sub_bgn_pos = sub_bgn_pos + N/2;
    sub_end_pos = sub_end_pos + N/2;
  }
  while (!last);
  return fasta_seq_vec;
}

}  // namespace fasta_sub_util

}  // namespace toppic
