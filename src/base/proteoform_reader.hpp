#ifndef PROT_BASE_PROTEOFORM_READER_HPP_
#define PROT_BASE_PROTEOFORM_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

#include "htslib/faidx.h"

#include "base/string_util.hpp"
#include "base/residue_seq.hpp"
#include "base/fasta_reader.hpp"
#include "base/proteoform.hpp"

namespace prot {

class ProteoformReader {
 public:
  /**
   * Constructs an instance with a File.
   **/
  ProteoformReader(const std::string &file_name);

  ProteoformPtr getNextProteoformPtr(const ResiduePtrVec &residue_list);

  void setSeqId (int id) {seq_id_ = id;}

  static ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                                const ResiduePtrVec &residue_list);


  static ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                                const ResiduePtrVec &residue_list,
                                                int seq_bgn_id);

 private:
  FastaReaderPtr reader_ptr_;
  int seq_id_ = 0;
};

}
#endif
