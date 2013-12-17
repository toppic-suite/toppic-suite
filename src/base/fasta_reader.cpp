#include <base/logger.hpp>

#include "base/fasta_reader.hpp"

namespace prot {

FastaReader::FastaReader(std::string file_name) {
  input_.open(file_name.c_str(), std::ios::in);
  if (!input_.is_open()) {
    LOG_ERROR( "fasta file  " << file_name << " does not exist.");
    throw "fasta file does not exist.";
  }
  std::getline(input_, ori_name_);
}

std::vector<std::string> FastaReader::getNextSeq() {
  std::vector<std::string> strs;
  if (!input_.is_open()) {
    return strs;
  }

  /* get the letters of sequence */
  std::string ori_seq;
  std::string prot_name = ori_name_.substr(1, ori_name_.size() -1 );
  std::string line;
  while (std::getline(input_, line)) {
    if (line.length() >= 1 && line.substr(0, 1) == ">") {
      ori_name_ = line;
      return fastaPreprocess(prot_name, ori_seq);
    }
    line = trim(line);
    ori_seq = ori_seq + line;
  }
  input_.close();
  return fastaPreprocess(prot_name, ori_seq);
}

/**
 * Read FASTA file and return next protein as an ResSeq.
 **/
ProteoformPtr FastaReader::getNextProteoformPtr(AcidPtrVec acid_list, 
                                             ResiduePtrVec residue_list) {
  std::vector<std::string> seq_info = getNextSeq();
  if (seq_info.size() == 0) {
    return ProteoformPtr(nullptr);
  }
  std::string name = seq_info[0];
  std::string seq = seq_info[1];
  LOG_TRACE( "name " << seq_info[0] << " seq " << seq_info[1]);
  AcidPtrVec acid_seq = convertSeqToAcidSeq(acid_list, seq); 
  ResiduePtrVec residue_ptrs = convertAcidToResidueSeq(residue_list, acid_seq);
  ResSeqPtr residue_seq_ptr = ResSeqPtr(new ResidueSeq(residue_ptrs)); 
  return getOriProteoformPtr(name, residue_seq_ptr);
}

/** process fasta string and remove unknown letters */
std::string rmChar(std::string str) {
  std::string seq = "";
  for (unsigned int i = 0; i < str.length(); i++) {
    char c = str.at(i);
    if (c < 'A' || c > 'Z') {
      continue;
    }
    char r = c;
    if (c == 'B') {
      r = 'D';
    } else if (c == 'Z') {
      r = 'E';
    } else if (c == 'X') {
      r = 'A';
    }
    seq = seq + r;
  }
  return seq;
}

/** Process the string */
std::vector<std::string> fastaPreprocess(std::string name, std::string seq) {
  std::string new_seq = rmChar(seq);
  if (!(new_seq == seq)) {
    LOG_INFO( "Reading sequence. Unknown letter occurred. ");
  }
  std::vector<std::string> strs;
  strs.push_back(name);
  strs.push_back(new_seq);
  return strs;
}

ProteoformPtrVec readFastaToProteoform(std::string file_name, 
                                       AcidPtrVec &acid_list, ResiduePtrVec &residue_list) {

  LOG_DEBUG( "start open file " << file_name);
  FastaReader reader(file_name);
  LOG_DEBUG( "open file done " << file_name);
  ProteoformPtrVec list;
  ProteoformPtr ptr;
  while ((ptr = reader.getNextProteoformPtr(acid_list, residue_list)).get() != nullptr) {
    list.push_back(ptr);
  }
  return list;
}

}

