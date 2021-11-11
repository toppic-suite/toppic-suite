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

#include "seq/fasta_reader.hpp"
#include "search/duplicatematch/additional_match.hpp"
#include "prsm/prsm.hpp"

namespace toppic {
  std::string AdditionalMatch::trimSequence(std::string raw_prot) {
    std::string new_seq = raw_prot;
    std::size_t found_1 = new_seq.find("(");
    std::size_t found_2 = new_seq.find(")");
    
    while (found_1 != std::string::npos || found_2 != std::string::npos) {// if ")" or "(" in the seq
        std::string left_seq = new_seq.substr(0, found_1);
        std::string right_seq = new_seq.substr(found_2 + 1);

        std::string in_parenthesis = new_seq.substr(found_1 + 1, found_2 - found_1 -1);

        new_seq = left_seq + in_parenthesis + right_seq;

        found_1 = new_seq.find("(");
        found_2 = new_seq.find(")");
    }

    found_1 = new_seq.find("[");
    found_2 = new_seq.find("]");

    while (found_1 != std::string::npos || found_2 != std::string::npos) {
        std::string left_seq = new_seq.substr(0, found_1);
        std::string right_seq = new_seq.substr(found_2 + 1);
        
        new_seq = left_seq + right_seq;

        found_1 = new_seq.find("[");
        found_2 = new_seq.find("]");
    }
    found_1 = new_seq.find(".");

    if (found_1 != std::string::npos) {
      new_seq = new_seq.substr(found_1 + 1);
    }
    found_1 = new_seq.find(".");
    
    if (found_1 != std::string::npos) {
      new_seq = new_seq.substr(0, found_1);
    }   
    return new_seq; 
  }

  void AdditionalMatch::process(PrsmPtr prsm_ptr_) {
    FastaReader reader(db_file_name_);
    FastaSeqPtr seq_ptr = reader.getNextSeq();
    std::vector<PrsmPtr> additional_match_prsms;

    std::string prot_seq = trimSequence(prsm_ptr_->getProteoformPtr()->getProteinMatchSeq());
    bool isExactMatchFound = false;
    int hit_cnt = 0;

    while (seq_ptr != nullptr) {
      std::string fasta_seq = seq_ptr->getRawSeq();
      std::size_t found = fasta_seq.find(prot_seq);
      if (found != std::string::npos) {
        isExactMatchFound = true;
        hit_cnt++;

        PrsmPtr additional_match_prsm_ptr = std::make_shared<Prsm>(*prsm_ptr_);
        additional_match_prsm_ptr->setPrsmId(prsm_cnt_);
        additional_match_prsm_ptr->setIsExactMatch(true);
        additional_match_prsm_ptr->setHitCnt(hit_cnt); // +1 is for the original match
        additional_match_prsm_ptr->getProteoformPtr()->setStartPos(found + 1);
        additional_match_prsm_ptr->getProteoformPtr()->setEndPos(found + prot_seq.size());
        additional_match_prsm_ptr->getProteoformPtr()->setFastaSeqPtr(seq_ptr);

        additional_match_prsms.push_back(additional_match_prsm_ptr);
        prsm_cnt_++;
      }
      seq_ptr = reader.getNextSeq();
    }
    if (!isExactMatchFound) {
      //search for approximate matches
      std::cout << prsm_ptr_->getProteoformPtr()->getProteinMatchSeq() << std::endl;
    }

    prsm_ptr_vec_.insert(prsm_ptr_vec_.end(), additional_match_prsms.begin(), additional_match_prsms.end());
  }
}