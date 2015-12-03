#include "base/logger.hpp"
#include "base/acid_util.hpp"
#include "base/residue_util.hpp"
#include "base/proteoform_reader.hpp"
#include "base/proteoform_factory.hpp"

namespace prot {

/*
ProteoformPtr getNextProteoform (FastaSeqPtr seq_ptr, FixModPtrVec fix_mod_list) {
  LOG_TRACE("name " << seq_ptr->getName() << " seq " << seq_ptr->getSeq());
  ResiduePtrVec residue_ptrs = ResidueUtil::convertSeqToResiduePtrVec(seq_ptr->getSeq());
  
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    for (size_t j = 0; j < fix_mod_list.size(); j++) {
      if (residue_ptrs[i] == fix_mod_list[j]->getOriResiduePtr()) {
        residue_ptrs[i] = fix_mod_list[j]->getModResiduePtr();
        break;
    }
  }
  
  DbResSeqPtr db_residue_seq_ptr(
      new DbResidueSeq(residue_ptrs, seq_ptr->getName(), seq_ptr->getDesc())); 
  seq_id_++;
  return ProteoformFactory::geneDbProteoformPtr(db_residue_seq_ptr);

ProteoformPtrVec ProteoformReader::readFastaToProteoform(const std::string &file_name,
                                                         const FixModPtrVec &fix_mod_list) {
  return readFastaToProteoform(file_name, residue_list, 0);
}


ProteoformPtrVec ProteoformReader::readFastaToProteoform(const std::string &file_name, 
                                                         const FixModPtrVec &fix_mod_list,
                                                         int seq_bgn_id) {
  LOG_DEBUG( "start open file " << file_name);
  FastaSeqReader reader(file_name);
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
*/

}

