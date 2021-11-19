//Copyright (c) 2014 - 2021, The Trustees of Indiana University.
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

#include <iostream>
#include <string>  
#include <vector>  
#include <sstream>

#include "common/util/logger.hpp"
#include "seq/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/search_fasta_match.hpp"

namespace toppic {

SearchFastaMatch::SearchFastaMatch(std::string db_file_name):
  db_file_name_(db_file_name) {
    FastaReaderPtr reader_ptr = std::make_shared<FastaReader>(db_file_name_);
    FastaSeqPtr seq_ptr = reader_ptr->getNextSeq();
    while (seq_ptr != nullptr) {
      fasta_seq_vec_.push_back(seq_ptr);
      seq_vec_.push_back(seq_ptr->getRawSeq());
      seq_ptr = reader_ptr->getNextSeq();
    }
  }

std::vector<std::pair<FastaSeqPtr, int>> SearchFastaMatch::process(PrsmPtr prsm_ptr_) {
  std::string raw_seq = prsm_ptr_->getProteoformPtr()->getFastaSeqPtr()->getRawSeq();
  int start_pos = prsm_ptr_->getProteoformPtr()->getStartPos();
  int end_pos = prsm_ptr_->getProteoformPtr()->getEndPos();

  std::string form_seq = raw_seq.substr(start_pos, end_pos - start_pos + 1);
  LOG_DEBUG("proteoform seq " << form_seq);

  std::vector<std::pair<FastaSeqPtr, int>> matches;
  for (size_t i = 0; i < seq_vec_.size(); i++) {
    std::size_t found = seq_vec_[i].find(form_seq);
    if (found != std::string::npos) {
      std::pair<FastaSeqPtr, int> match(fasta_seq_vec_[i], found);
      matches.push_back(match);
    }
  }
  return matches;
}

}
