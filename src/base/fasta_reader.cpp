
#include "base/logger.hpp"
#include "base/fasta_reader.hpp"
#include "base/db_block.hpp"

namespace prot {

FastaSeq::FastaSeq(const std::string &name_line, const std::string &ori_seq) {
  int space_pos = name_line.find(" ");
  name_ = name_line.substr(0, space_pos);
  desc_ = name_line.substr(space_pos + 1);
  seq_ = rmChar(ori_seq);
  //LOG_DEBUG("name line:" << name_line);
  //LOG_DEBUG("NAME:" << name_ << "DESC:" << desc_);
}

/** process fasta string and remove unknown letters */
std::string rmChar(const std::string &ori_seq) {
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
      new DbResidueSeq(residue_ptrs, seq_id_, seq_ptr->getName(), seq_ptr->getDesc())); 
  seq_id_++;
  return getDbProteoformPtr(db_residue_seq_ptr);
}

ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list) {
  return readFastaToProteoform(file_name, residue_list, 0);
}

ProteoformPtrVec readFastaToProteoform(const std::string &file_name, 
                                       const ResiduePtrVec &residue_list,
                                       int seq_bgn_id) {
  LOG_DEBUG( "start open file " << file_name);
  FastaReader reader(file_name);
  reader.setSeqId (seq_bgn_id);
  
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

ProteoformPtr readFastaToProteoform(faidx_t *fai,
                                    int id,
                                    const std::string &seq_name, 
                                    const std::string &seq_desc,
                                    const ResiduePtrVec &residue_list) {
  int seq_len;
  char *seq = fai_fetch(fai, seq_name.c_str(), &seq_len);
  if ( seq_len < 0 ) {
    LOG_ERROR("Failed to fetch sequence " << seq_name);
  }
  std::string raw_seq_str(seq);
  free(seq);
  std::string seq_str = rmChar(raw_seq_str);
  //LOG_DEBUG("Seq name " << seq_name << " Seq: " << seq_str);
  AcidPtrVec acid_seq = AcidFactory::convertSeqToAcidSeq(seq_str); 
  ResiduePtrVec residue_ptrs = convertAcidToResidueSeq(residue_list, acid_seq);
  DbResSeqPtr db_residue_seq_ptr(
      new DbResidueSeq(residue_ptrs, id, seq_name, seq_desc)); 
  return getDbProteoformPtr(db_residue_seq_ptr);
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
    std::string desc = seq_info->getDesc();
    std::string decoy_name = "DECOY_" + name;
    std::string temp = seq.substr(2, seq.length() - 2);
    std::random_shuffle(temp.begin(), temp.end());
    std::string decoy_seq = seq.substr(0,2) + temp;
    output << ">" << decoy_name << " " << desc <<  std::endl;
    output << decoy_seq << std::endl;
    output << ">" << name << " " << desc << std::endl;
    output << seq << std::endl;

    seq_info = reader.getNextSeq();
  }
  output.close();
}

void generateStandardDb(const std::string &ori_file_name, 
                        const std::string &st_file_name) {

  std::ifstream ori_db(ori_file_name);
  std::string line;
  std::ofstream standard_db;
  standard_db.open(st_file_name.c_str(), std::ios::out);

  while (std::getline(ori_db, line)) {
    if(line.length() > 0){
      standard_db << line << std::endl;
    }  
  }
  ori_db.close();
  standard_db.close();
}

void generateDbBlock(const std::string &db_file_name, 
                     int block_size) {
  int block_idx = 0;
  int seq_idx = 0;

  std::ofstream index_output;
  std::string index_file_name = db_file_name + "_block_index";
  index_output.open(index_file_name.c_str(), std::ios::out);
  std::ofstream block_output;
  std::string block_file_name = db_file_name + "_" + std::to_string(block_idx);
  block_output.open(block_file_name.c_str(), std::ios::out);

  FastaReader reader(db_file_name);
  FastaSeqPtr seq_info = reader.getNextSeq();
  index_output << block_idx << "\t" << seq_idx << std::endl;
  int seq_size = 0;
  while (seq_info!=nullptr) {
    std::string name = seq_info->getName();
    std::string desc = seq_info->getDesc();
    std::string seq = seq_info->getSeq();
    seq_size += seq.length();
    if (seq_size > block_size) {
      block_output.close();
      block_idx++;
      LOG_DEBUG("Database block " << block_idx << " size " << seq_size);
      index_output << block_idx << "\t" << seq_idx << std::endl;
      seq_size = 0;
      block_file_name = db_file_name + "_" + std::to_string(block_idx);
      block_output.open(block_file_name.c_str(), std::ios::out);
    }
    block_output << ">" << name << " " << desc << std::endl;
    block_output << seq << std::endl;
    seq_info = reader.getNextSeq();
    seq_idx++;
  }
  index_output.close();
  block_output.close();
}
  

void dbPreprocess(const std::string &ori_db_file_name, 
                  const std::string &db_file_name, 
                  bool decoy, int block_size) {

  std::string standard_db_file_name = ori_db_file_name + "_standard";
  generateStandardDb(ori_db_file_name, standard_db_file_name);
  
  if (decoy) {
    generateShuffleDb(standard_db_file_name, db_file_name);
  }
  else {
    boost::filesystem::path ori_path(standard_db_file_name);
    boost::filesystem::path db_path(db_file_name);
    boost::filesystem::copy_file(ori_path, db_path, 
                                 boost::filesystem::copy_option::overwrite_if_exists);
  }

  generateDbBlock(db_file_name, block_size);

  fai_build(db_file_name.c_str());
}

}

