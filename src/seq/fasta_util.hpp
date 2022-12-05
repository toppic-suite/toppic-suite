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

#ifndef TOPPIC_SEQ_FASTA_UTIL_HPP_
#define TOPPIC_SEQ_FASTA_UTIL_HPP_

#include <string>

#include "common/util/str_util.hpp"

namespace toppic {

namespace fasta_util {

std::string getString(const StringPairVec &str_pair_vec);

void generateShuffleDb(const std::string &file_name,
                       const std::string &target_decoy_file_name);

void dbSimplePreprocess(const std::string &ori_db_file_name,
                        const std::string &db_file_name);

void dbPreprocess(const std::string &ori_db_file_name,
                  const std::string &db_file_name,
                  bool decoy, int block_size, 
                  int max_frag_len, int min_block_num);

int countProteinNum(const std::string &fasta_file);

}  // namespace fasta_util

}  // namespace toppic

#endif
