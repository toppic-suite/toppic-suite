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
#include "prsm/prsm_para.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class TableWriter {
 public:
  TableWriter(PrsmParaPtr prsm_para_ptr, 
              const std::string &input_file_ext, 
              const std::string &output_file_ext);
  void write();
 private:
  PrsmParaPtr prsm_para_ptr_;
  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<TableWriter> TableWriterPtr;

} /* namespace prot */

#endif /* TABLE_WRITER_HPP_ */
