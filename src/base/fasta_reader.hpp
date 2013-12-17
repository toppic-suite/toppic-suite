#ifndef PROT_FASTA_READER_HPP_
#define PROT_FASTA_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

#include "base/residue_seq.hpp"
#include "base/proteoform.hpp"

namespace prot {

class FastaReader {
 public:
  /**
   * Constructs an instance with a File.
   **/
  FastaReader(std::string);

  /**
   * Read FASTA file and return next protein
   * name and sequence. result[0] is
   * protein name and result[1] is sequence.
   **/                  
  std::vector<std::string> getNextSeq();

  ProteoformPtr getNextProteoformPtr(AcidPtrVec acid_list, 
                                  ResiduePtrVec residue_list);
 private:
  std::ifstream input_;
  std::string ori_name_;
};

std::vector<std::string> fastaPreprocess(std::string name, std::string seq);

ProteoformPtrVec readFastaToProteoform(std::string file_name, 
                                       AcidPtrVec &acid_list, ResiduePtrVec &residue_list);

inline std::string &trim(std::string &s) {
  s.erase(s.begin(), 
          std::find_if(s.begin(), s.end(), 
                       std::not1(std::ptr_fun<int, int>(std::isspace))));
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))).base(), 
          s.end());
  return s;
}

}

#endif
