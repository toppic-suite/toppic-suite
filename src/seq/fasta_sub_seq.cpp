// Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

namespace toppic {

FastaSubSeq::FastaSubSeq(const std::string &name_line, const std::string &sub_seq,
                         int sub_seq_start):
    FastaSeq(name_line, sub_seq),
    sub_seq_start_(sub_seq_start) {
      length_ = static_cast<int>(sub_seq.length());
      sub_seq_end_ = sub_seq_start_ + length_ - 1;
    }

FastaSubSeq::FastaSubSeq(const std::string &name, const std::string &desc,
                         const std::string &sub_seq, int sub_seq_start):
    FastaSeq(name, desc, sub_seq),
    sub_seq_start_(sub_seq_start) {
      length_ = static_cast<int>(sub_seq.length());
      sub_seq_end_ = sub_seq_start_ + length_ - 1;
    }

}  // namespace toppic
