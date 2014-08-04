#include "base/logger.hpp"
#include "base/fasta_reader.hpp"

namespace prot {

/** process fasta string and remove unknown letters */
std::string FastaSeq::rmChar(const std::string &ori_seq) {
  std::string seq = "";
  for (size_t i = 0; i < ori_seq.length(); i++) {
    char c = ori_seq.at(i);
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
  if (ori_seq != seq) { 
    LOG_INFO( "Reading sequence. Unknown letter occurred. ");
  }
  return seq;
}

FastaReader::FastaReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
  if (!input_.is_open()) {
    LOG_ERROR( "fasta file  " << file_name << " does not exist.");
    throw "fasta file does not exist.";
  }
  std::getline(input_, ori_name_);
}

FastaSeqPtr FastaReader::getNextSeq() {
  if (!input_.is_open()) {
    return FastaSeqPtr(nullptr);
  }

  /* get the letters of sequence */
  std::string ori_seq;
  std::string prot_name = ori_name_.substr(1, ori_name_.size() -1 );
  std::string line;
  while (std::getline(input_, line)) {
    if (line.length() >= 1 && line.substr(0, 1) == ">") {
      ori_name_ = trim(line);
      return FastaSeqPtr(new FastaSeq(prot_name, ori_seq));
    }
    line = trim(line);
    ori_seq = ori_seq + line;
    if (ori_seq.size() >= 1000000) {
      LOG_ERROR("Protein sequences are too long! Incorrect fasta file!");
      throw("fasta file error");
    }
  }
  input_.close();
  return FastaSeqPtr(new FastaSeq(prot_name, ori_seq));
}
/**
 * Read FASTA file and return next protein as an ResSeq.
 * residue_list determine fixed PTMs
 **/
ProteoformPtr FastaReader::getNextProteoformPtr(const ResiduePtrVec &residue_list) {
  FastaSeqPtr seq_ptr = getNextSeq();
  if (seq_ptr.get() == nullptr) {
    return ProteoformPtr(nullptr);
  }
  LOG_TRACE("name " << seq_ptr->getName() << " seq " << seq_ptr->getSeq());
  AcidPtrVec acid_seq = AcidFactory::convertSeqToAcidSeq(seq_ptr->getSeq()); 
  ResiduePtrVec residue_ptrs = convertAcidToResidueSeq(residue_list, acid_seq);
  DbResSeqPtr db_residue_seq_ptr(
      new DbResidueSeq(residue_ptrs, id_, seq_ptr->getName())); 
  id_++;
  return getDbProteoformPtr(db_residue_seq_ptr);
}

ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list) {
  LOG_DEBUG( "start open file " << file_name);
  FastaReader reader(file_name);
  LOG_DEBUG( "open file done " << file_name);
  ProteoformPtrVec list;
  ProteoformPtr ptr = reader.getNextProteoformPtr(residue_list);
  int count = 0;
  while (ptr.get() != nullptr) {
    list.push_back(ptr);
    ptr = reader.getNextProteoformPtr(residue_list);
    count++;
  }
  return list;
}

/** initialize sequence list */
void generateShuffleDb(const std::string &file_name, 
                       const std::string &target_decoy_file_name) {

  std::ofstream output;
  output.open(target_decoy_file_name.c_str(), std::ios::out);
  FastaReader reader(file_name);

  FastaSeqPtr seq_info = reader.getNextSeq();
  while (seq_info!=nullptr) {
    std::string name = seq_info->getName();
    std::string seq = seq_info->getSeq();
    std::string decoy_name = "DECOY_" + name;
    std::string temp = seq.substr(2, seq.length() - 2);
    std::random_shuffle(temp.begin(), temp.end());
    std::string decoy_seq = seq.substr(0,2) + temp;
    output << ">" << decoy_name << std::endl;
    output << decoy_seq << std::endl;
    output << ">" << name << std::endl;
    output << seq << std::endl;

    seq_info = reader.getNextSeq();
  }
  output.close();
}

}

