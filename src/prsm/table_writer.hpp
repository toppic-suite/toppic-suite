/*
 * table_writer.hpp
 *
 *  Created on: Feb 19, 2014
 *      Author: xunlikun
 */

#ifndef TABLE_WRITER_HPP_
#define TABLE_WRITER_HPP_

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <map>

#include "base/string_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/species.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class TableWriter {
 public:
  TableWriter(std::string db_file,std::string spec_file,std::string input_file,std::string output_file,double ppo);
  TableWriter(std::map<std::string,std::string> arguments,
                           std::string input_file,
                           std::string output_file);
  void write();
 private:
  std::string spec_file_;
  std::string db_file_;
  std::string input_file_;
  std::string output_file_;
  std::ofstream file_;
  double ppo_;
};

typedef std::shared_ptr<TableWriter> TableWriterPtr;

} /* namespace prot */

#endif /* TABLE_WRITER_HPP_ */
