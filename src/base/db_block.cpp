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

DbBlockPtrVec readDbBlockIndex(const std::string &db_file_name) {
  DbBlockPtrVec db_block_ptr_vec;
  std::ifstream input;
  std::string index_file_name = db_file_name + "_index";
  input.open(index_file_name.c_str(), std::ios::in);
  if (!input.is_open()) {
    LOG_ERROR( "index file  " << index_file_name << " does not exist.");
    throw "index file does not exist.";
  }
  std::string line;
  std::vector<std::string> strs;
  while (std::getline(input, line)) {
      strs = split(line, '\t');
      int block_index = std::stoi(strs[0]);
      int seq_index = std::stoi(strs[1]);
      LOG_DEBUG("block " << block_index << " seq " << seq_index);
      DbBlockPtr ptr(new DbBlock(block_index, seq_index));
      db_block_ptr_vec.push_back(ptr);
  }
  input.close();
  return db_block_ptr_vec;
}

void generateDbBlockIndex(const std::string &db_file_name, 
                          int block_size) {
  int block_idx = 0;
  int seq_idx = 0;

  std::ofstream index_output;
  std::string index_file_name = db_file_name + "_block_index";
  index_output.open(index_file_name.c_str(), std::ios::out);
  std::ofstream block_output;
  std::string block_file_name = db_file_name + "_" + std::to_string(block_idx);
  block_output.open(block_file_name.c_str(), std::ios::out);

  FastaReader reader(db_file_name);
  FastaSeqPtr seq_info = reader.getNextSeq();
  index_output << block_idx << "\t" << seq_idx << std::endl;
  int seq_size = 0;
  while (seq_info!=nullptr) {
    std::string name = seq_info->getName();
    std::string seq = seq_info->getSeq();
    seq_size += seq.length();
    if (seq_size > block_size) {
      block_output.close();
      block_idx++;
      LOG_DEBUG("Database block " << block_idx << " size " << seq_size);
      index_output << block_idx << "\t" << seq_idx << std::endl;
      seq_size = 0;
      block_file_name = db_file_name + "_" + std::to_string(block_idx);
      block_output.open(block_file_name.c_str(), std::ios::out);
    }
    block_output << ">" << name << std::endl;
    block_output << seq << std::endl;
    seq_info = reader.getNextSeq();
    seq_idx++;
  }
  index_output.close();
  block_output.close();
}

}
