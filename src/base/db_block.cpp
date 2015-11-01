#include "base/db_block.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "base/fasta_reader.hpp"

namespace prot {

DbBlock::DbBlock(int block_index, int seq_index) {
  block_index_ = block_index;
  seq_index_ = seq_index;
}

DbBlockPtrVec DbBlock::readDbBlockIndex(const std::string &db_file_name) {
  DbBlockPtrVec db_block_ptr_vec;
  std::ifstream input;
  std::string index_file_name = db_file_name + "_block_index";
  input.open(index_file_name.c_str(), std::ios::in);
  if (!input.is_open()) {
    LOG_ERROR( "index file  " << index_file_name << " does not exist.");
    throw "index file does not exist.";
  }
  std::string line;
  std::vector<std::string> strs;
  while (std::getline(input, line)) {
    strs = StringUtil::split(line, '\t');
    int block_index = std::stoi(strs[0]);
    int seq_index = std::stoi(strs[1]);
    LOG_DEBUG("block " << block_index << " seq " << seq_index);
    DbBlockPtr ptr(new DbBlock(block_index, seq_index));
    db_block_ptr_vec.push_back(ptr);
  }
  input.close();
  return db_block_ptr_vec;
}

}
