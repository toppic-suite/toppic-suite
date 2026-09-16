// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "seq/fasta_index_reader.hpp"

#include <cstdlib>
#include <memory>
#include <mutex>
#include <string>

#include "common/util/logger.hpp"

namespace toppic {

namespace {
// fai_fetch is not thread safe, so all fetches are serialized on this mutex.
std::mutex fasta_index_reader_mtx;
}  // namespace

FastaIndexReader::FastaIndexReader(const std::string& file_name) {
  fai_ = fai_load(file_name.c_str());
}

FastaIndexReader::~FastaIndexReader() { fai_destroy(fai_); }

FastaSeqPtr FastaIndexReader::readFastaSeq(const std::string& name,
                                           const std::string& desc) {
  std::string ori_seq;
  {
    std::lock_guard<std::mutex> lock(fasta_index_reader_mtx);
    int seq_len = 0;
    char* seq = fai_fetch(fai_, name.c_str(), &seq_len);
    if (seq == nullptr || seq_len < 0) {
      LOG_WARN("Failed to fetch protein sequence " << name);
    } else {
      ori_seq.assign(seq, seq_len);
    }
    free(seq);
  }
  return std::make_shared<FastaSeq>(name, desc, ori_seq);
}

}  // namespace toppic
