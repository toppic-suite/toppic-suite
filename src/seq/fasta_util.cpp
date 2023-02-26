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

#include <string>
#include <algorithm>
#include <random>
#include <cmath>

#include "htslib/faidx.h"

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/base/ptm_base.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/fasta_util.hpp"

namespace toppic {

namespace fasta_util {

std::string getString(const StringPair &str_pair) {
  std::string result = str_pair.first;
  std::string ptm_str = str_pair.second;
  if (ptm_str != PtmBase::getEmptyPtmPtr()->getAbbrName()) {
    result = result + "[" + ptm_str + "]";
  }
  return result;
}

std::string getString(const StringPairVec &str_pair_vec) {
  std::string result;
  for (size_t i = 0; i < str_pair_vec.size(); i++) {
    result = result + getString(str_pair_vec[i]);
  }
  return result;
}

void generateShuffleDb(const std::string &file_name,
                       const std::string &target_decoy_file_name) {
  std::ofstream output;
  output.open(target_decoy_file_name.c_str(), std::ios::out);
  FastaReader reader(file_name);
  FastaSeqPtr seq_info = reader.getNextSeq();
  std::mt19937 r{std::random_device{}()};
  r.seed(std::mt19937::default_seed);
  while (seq_info != nullptr) {
    std::string name = seq_info->getName();
    std::string seq = seq_info->getRawSeq();
    StringPairVec str_pair_vec = seq_info->getAcidPtmPairVec();
    std::string desc = seq_info->getDesc();
    std::string decoy_name = "DECOY_" + name;
    std::string decoy_seq;
    if (str_pair_vec.size() > 2) {
      std::shuffle(str_pair_vec.begin() + 2, str_pair_vec.end(), r);
      decoy_seq = getString(str_pair_vec);
    } else {
      decoy_seq = seq;
    }
    output << ">" << decoy_name << " " << desc <<  std::endl;
    output << decoy_seq << std::endl;
    output << ">" << name << " " << desc << std::endl;
    output << seq << std::endl;

    seq_info = reader.getNextSeq();
  }
  output.close();
}

void generateStandardDb(const std::string &ori_file_name,
                        const std::string &st_file_name) {
  std::ifstream ori_db(ori_file_name);
  std::string line;
  std::ofstream standard_db;

  standard_db.open(st_file_name.c_str(), std::ios::out);

  while (std::getline(ori_db, line)) {
    if (line.length() > 0) {
      standard_db << line << std::endl;
    }
  }
  ori_db.close();
  standard_db.close();
}

void generateDbBlock(const std::string &db_file_name, int block_size, 
                     int max_frag_len, int min_block_num) {
  int block_idx = 0;
  int seq_idx = 0;
  std::string index_file_name = db_file_name + "_block_index";
  std::string block_file_name = db_file_name + "_" + str_util::toString(block_idx);


  // get total_size;
  FastaReader size_reader(db_file_name);
  FastaSeqPtr seq_info = size_reader.getNextSeq();
  double total_seq_size = 0;
  int seq_num = 0;
  while (seq_info != nullptr) {
    std::string seq = seq_info->getRawSeq();
    int seq_len = seq.length();
    if (seq_len <= max_frag_len) {
      total_seq_size = total_seq_size + (seq_len * (seq_len + 1)/2);
    }
    else {
      total_seq_size = total_seq_size + (seq_len - max_frag_len) * max_frag_len 
	      + (max_frag_len * (max_frag_len + 1) /2);
    }
    seq_num++;
    seq_info = size_reader.getNextSeq();
  }
  size_reader.close();
  //LOG_ERROR("total seq size " << total_seq_size); 
  // adjust block size if the database is too small and the sequence number 
  // is large enough > 50 * min_block_num
  int block_num = static_cast<int>(std::ceil(total_seq_size/block_size));
  if (block_num < min_block_num && seq_num > (min_block_num * 50)) {
    block_num = min_block_num;
    // add +1 to make sure there is no error introduced in division.
    block_size = total_seq_size/block_num + 1;
  }

  std::ofstream index_output;
  index_output.open(index_file_name.c_str(), std::ios::out);
  std::ofstream block_output;
  block_output.open(block_file_name.c_str(), std::ios::out);

  FastaReader reader(db_file_name);
  seq_info = reader.getNextSeq();
  index_output << block_idx << "\t" << seq_idx << std::endl;
  long seq_size = 0;
  while (seq_info != nullptr) {
    std::string name = seq_info->getName();
    std::string desc = seq_info->getDesc();
    std::string seq = seq_info->getRawSeq();
    int seq_len = seq.length();
    if (seq_len <= max_frag_len) {
      seq_size = seq_size + (seq_len * (seq_len + 1)/2);
    }
    else {
      seq_size = seq_size + (seq_len - max_frag_len) * max_frag_len 
	      + (max_frag_len * (max_frag_len + 1) /2);
    }
    if (seq_size > block_size) {
      block_output.close();
      block_idx++;
      LOG_DEBUG("Database block " << block_idx << " size " << seq_size);
      index_output << block_idx << "\t" << seq_idx << std::endl;
      seq_size = 0;
      block_file_name = db_file_name + "_" + str_util::toString(block_idx);
      block_output.open(block_file_name.c_str(), std::ios::out);
    }
    block_output << ">" << name << " " << desc << std::endl;
    block_output << seq << std::endl;
    seq_info = reader.getNextSeq();
    seq_idx++;
  }
  index_output.close();
  block_output.close();
}

void dbSimplePreprocess(const std::string &ori_db_file_name,
                        const std::string &db_file_name) {
  if (!file_util::exists(ori_db_file_name + "_idx")){//if _idx folder doesn't exist yet
    file_util::createFolder(ori_db_file_name + "_idx");
  }
  if (file_util::exists(db_file_name)) {
    return;
  }
  generateStandardDb(ori_db_file_name, db_file_name);

  if (file_util::exists(db_file_name + ".fai")) {
    return;
  }
  fai_build(db_file_name.c_str());
}

void dbPreprocess(const std::string &ori_db_file_name,
                  const std::string &file_name,
                  bool decoy, int block_size, 
                  int max_frag_len, int min_block_num) {
  //if _idx folder doesn't exist yet
  if (!file_util::exists(ori_db_file_name + "_idx")){
    file_util::createFolder(ori_db_file_name + "_idx");
  }
  // Generate a stardard fasta file in which empty lines are removed
  std::string standard_db_file_name = ori_db_file_name + "_idx" 
    + file_util::getFileSeparator() 
    + file_util::filenameFromEntirePath(ori_db_file_name) + "_standard";  

  // Generate new database file
  std::string new_db_file_name = ori_db_file_name + "_idx" 
    + file_util::getFileSeparator() 
    + file_util::filenameFromEntirePath(file_name);
  
  if (standard_db_file_name.size() > 200 || new_db_file_name.size() > 200) {
    LOG_ERROR("Database file name is too long. Please rename it and try running again.");
    exit(1);
  }
  if (!file_util::exists(standard_db_file_name)) {
    generateStandardDb(ori_db_file_name, standard_db_file_name);
  }

  if (decoy) {
    if (!file_util::exists(new_db_file_name)) {
      generateShuffleDb(standard_db_file_name, new_db_file_name);
    }
  } 
  else {
    if (!file_util::exists(new_db_file_name)) {
      bool overwrite = true;
      file_util::copyFile(standard_db_file_name, new_db_file_name, overwrite); 
    }
  }

  // Generate new database blocks
  if (!file_util::exists(new_db_file_name + "_block_index") || 
      !file_util::exists(new_db_file_name + "_0")) {
    generateDbBlock(new_db_file_name, block_size, max_frag_len, min_block_num);
  }

  // Generate new database file index
  if (!file_util::exists(new_db_file_name + ".fai")) {
    fai_build(new_db_file_name.c_str());
  }
}

int countProteinNum(const std::string &fasta_file) {
  FastaReader reader(fasta_file);
  int cnt = 0;
  while (reader.getNextSeq() != nullptr) {
    cnt++;
  }
  reader.close();
  return cnt;
}

}  // namespace fasta_util

}  // namespace toppic
