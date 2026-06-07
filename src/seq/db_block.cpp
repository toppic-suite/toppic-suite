//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#include "seq/db_block.hpp"

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"

namespace toppic {

DbBlock::DbBlock(int block_index, int seq_index):
      block_index_(block_index),
      seq_index_(seq_index) {}

DbBlockPtrVec DbBlock::readDbBlockIndex(const std::string &db_file_name) {
  std::string index_file_name = db_file_name + "_block_index";
  std::ifstream input(index_file_name);
  if (!input.is_open()) {
    LOG_ERROR("Index file " << index_file_name << " does not exist.");
    exit(EXIT_FAILURE);
  }
  DbBlockPtrVec db_block_ptr_vec;
  std::string line;
  while (std::getline(input, line)) {
    std::vector<std::string> strs = str_util::split(line, "\t");
    int block_index = std::stoi(strs[0]);
    int seq_index = std::stoi(strs[1]);
    db_block_ptr_vec.push_back(std::make_shared<DbBlock>(block_index, seq_index));
  }
  return db_block_ptr_vec;
}

}  // namespace toppic
