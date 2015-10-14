#ifndef PROT_PROTEOFORM_READER_HPP_
#define PROT_PROTEOFORM_READER_HPP_

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

 private:
  FastaReaderPtr reader_ptr_;
  int seq_id_ = 0;
};

ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list);
                                       

ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list,
                                       int seq_bgn_id);

ProteoformPtr readFastaToProteoform(faidx_t *fai,
                                    int id,
                                    const std::string &seq_name, 
                                    const std::string &seq_desc,
                                    const ResiduePtrVec &residue_list);

}

#endif
