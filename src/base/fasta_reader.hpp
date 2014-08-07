#ifndef PROT_FASTA_READER_HPP_
#define PROT_FASTA_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "base/string_util.hpp"
#include "base/residue_seq.hpp"
#include "base/proteoform.hpp"

namespace prot {

class FastaSeq {
 public:
  FastaSeq(const std::string &name, const std::string &ori_seq);

  std::string getName() {return name_;}

  std::string getSeq() {return seq_;}

 private:
  std::string name_;
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

  /**
   * Read FASTA file and return next protein
   * name and sequence. 
   **/                  
  FastaSeqPtr getNextSeq();

  ProteoformPtr getNextProteoformPtr(const ResiduePtrVec &residue_list);

 private:
  std::ifstream input_;
  std::string ori_name_;
  int id_ = 0;
};

ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list);

void generateShuffleDb(const std::string &file_name, 
                       const std::string &target_decoy_file_name);

}

#endif
