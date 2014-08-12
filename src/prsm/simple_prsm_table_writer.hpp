/*
 * simple_table_writer.hpp
 *
 *  Created on: Aug 11, 2014
 *      Author: qkou
 */

#ifndef PROT_SIMPLE_TABLE_WRITER_HPP_
#define PROT_SIMPLE_TABLE_WRITER_HPP_

#include "prsm/prsm_para.hpp"
#include "prsm/simple_prsm.hpp"

namespace prot {

class SimplePrsmTableWriter {
 public:
  SimplePrsmTableWriter(PrsmParaPtr prsm_para_ptr, 
                        const std::string &input_file_ext, 
                        const std::string &output_file_ext);
  void write();
 private:
  PrsmParaPtr prsm_para_ptr_;
  std::string input_file_ext_;
  std::string output_file_ext_;

};

typedef std::shared_ptr<SimplePrsmTableWriter> SimplePrsmTableWriterPtr;

}

#endif /* SIMPLE_TABLE_WRITER_HPP_ */
