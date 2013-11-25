#include <boost/algorithm/string.hpp>
#include "fasta_reader.hpp"

namespace prot {

FastaReader::FastaReader(const char* file_name) {
  input_.open(file_name, std::ios::in);
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
    std::cerr << "Reading sequence. Unknown letter occurred. " << seq;
  }
  std::vector<std::string> strs;
  strs.push_back(name);
  strs.push_back(new_seq);
  return strs;
}


}

