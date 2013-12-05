#include <boost/algorithm/string.hpp>

#include <log4cxx/logger.h>

#include "base/fasta_reader.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("FastaReader"));

FastaReader::FastaReader(std::string file_name) {
  input_.open(file_name.c_str(), std::ios::in);
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
    boost::algorithm::trim(line);
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
  LOG4CXX_TRACE(logger, "name " << seq_info[0] << " seq " << seq_info[1]);
  AcidPtrVec acid_seq = convertSeqToAcidSeq(acid_list, seq); 
  ResiduePtrVec residue_ptrs = convertAcidToResidueSeq(residue_list, acid_seq);
  ResSeqPtr residue_seq_ptr = ResSeqPtr(new ResidueSeq(residue_ptrs)); 
  return ProteoformPtr(new Proteoform(name, residue_seq_ptr));
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
    LOG4CXX_INFO(logger, "Reading sequence. Unknown letter occurred. ");
  }
  std::vector<std::string> strs;
  strs.push_back(name);
  strs.push_back(new_seq);
  return strs;
}

ProteoformPtrVec readFastaToProteoform(std::string file_name, 
                                       AcidPtrVec &acid_list, ResiduePtrVec &residue_list) {

  LOG4CXX_DEBUG(logger, "start open file " << file_name);
  FastaReader reader(file_name);
  LOG4CXX_DEBUG(logger, "open file done " << file_name);
  ProteoformPtrVec list;
  ProteoformPtr ptr;
  while ((ptr = reader.getNextProteoformPtr(acid_list, residue_list)).get() != nullptr) {
    list.push_back(ptr);
  }
  return list;
}

}

