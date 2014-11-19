#ifndef PROT_FASTA_READER_HPP_
#define PROT_FASTA_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS


#include "base/string_util.hpp"
#include "base/residue_seq.hpp"
#include "base/proteoform.hpp"

namespace prot {

class FastaSeq {
 public:
  FastaSeq(const std::string &name_line, const std::string &ori_seq);

  std::string getName() {return name_;}

  std::string getDesc() {return desc_;}

  std::string getSeq() {return seq_;}

 private:
  std::string name_;
  std::string desc_;
  std::string seq_;
  /* remove incorrect charaters in sequence */
  std::string rmChar(const std::string &ori_seq);
}; 

typedef std::shared_ptr<FastaSeq> FastaSeqPtr;

class FastaReader {
 public:
  /**
   * Constructs an instance with a File.
   **/
  FastaReader(const std::string &file_name);

  void setSeqId(int seq_id) {seq_id_ = seq_id;}

  /**
   * Read FASTA file and return next protein
   * name and sequence. 
   **/                  
  FastaSeqPtr getNextSeq();

  ProteoformPtr getNextProteoformPtr(const ResiduePtrVec &residue_list);

 private:
  std::ifstream input_;
  std::string ori_name_;
  int seq_id_ = 0;
};


ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list);
                                       

ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list,
                                       int seq_bgn_id);

void generateShuffleDb(const std::string &file_name, 
                       const std::string &target_decoy_file_name);

void dbPreprocess(const std::string &ori_db_file_name, 
                  const std::string &db_file_name, 
                  bool decoy, int block_size);

}

#endif
