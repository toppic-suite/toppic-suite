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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "seq/db_block.hpp"
#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "seq/fasta_reader.hpp"

namespace toppic {

DbBlockPtrVec DbBlock::readDbBlockIndex(const std::string &db_file_name) {
  DbBlockPtrVec db_block_ptr_vec;
  std::ifstream input;
  std::string index_file_name = db_file_name + "_block_index";
  input.open(index_file_name.c_str(), std::ios::in);
  if (!input.is_open()) {
    LOG_ERROR("index file  " << index_file_name << " does not exist.");
    throw "index file does not exist.";
  }
  std::string line;
  std::vector<std::string> strs;
  while (std::getline(input, line)) {
    strs = string_util::split(line, '\t');
    int block_index = std::stoi(strs[0]);
    int seq_index = std::stoi(strs[1]);
    LOG_DEBUG("block " << block_index << " seq " << seq_index);
    DbBlockPtr ptr = std::make_shared<DbBlock>(block_index, seq_index);
    db_block_ptr_vec.push_back(ptr);
  }
  input.close();
  return db_block_ptr_vec;
}

}  // namespace toppic
