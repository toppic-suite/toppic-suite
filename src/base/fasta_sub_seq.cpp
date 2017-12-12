// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <string>

#include "fasta_sub_seq.hpp"

namespace prot {

FastaSubSeq::FastaSubSeq(const std::string &name_line, const std::string &ori_seq,
                         int sub_seq_start):
    FastaSeq(name_line, ori_seq),
    sub_seq_start_(sub_seq_start) {
      length_ = static_cast<int>(ori_seq.length());
      sub_seq_end_ = sub_seq_start_ + length_ - 1;
    }

FastaSubSeq::FastaSubSeq(const std::string &name, const std::string &desc,
                         const std::string &ori_seq, int sub_seq_start):
    FastaSeq(name, desc, ori_seq),
    sub_seq_start_(sub_seq_start) {
      length_ = static_cast<int>(ori_seq.length());
      sub_seq_end_ = sub_seq_start_ + length_ - 1;
    }

}  // namespace prot
