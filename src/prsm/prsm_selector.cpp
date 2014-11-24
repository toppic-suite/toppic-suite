#include "base/file_util.hpp"
#include "prsm/prsm_selector.hpp"

namespace prot {

PrsmSelector::PrsmSelector(const std::string &db_file_name,
                           const std::string &spec_file_name,
                           const std::string &in_file_ext,
                           const std::string &out_file_ext, int n_top){
  spec_file_name_ = spec_file_name;
  db_file_name_ = db_file_name;
  input_file_ext_ = in_file_ext;
  output_file_ext_ = out_file_ext;
  n_top_ = n_top;
}

bool containsSameDbSeq(const PrsmPtrVec prsm_ptrs, PrsmPtr target_prsm_ptr) {
  for(size_t i=0; i< prsm_ptrs.size();i++){
    if(prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getId()==
        target_prsm_ptr->getProteoformPtr()->getDbResSeqPtr()->getId()){
      return true;
    }
  }
  return false;
}

PrsmPtrVec PrsmSelector::getTopPrsms(PrsmPtrVec &prsm_ptrs, int n_top){
  std::sort(prsm_ptrs.begin(),prsm_ptrs.end(),prsmEValueUp);
  int size = prsm_ptrs.size();
  int max = size > n_top? n_top:size;
  PrsmPtrVec result_ptrs;
  for(int i=0;i<max;i++){
    if(!containsSameDbSeq(result_ptrs, prsm_ptrs[i])){
      result_ptrs.push_back(prsm_ptrs[i]);
    }
  }
  return result_ptrs;
}

/*
void PrsmSelector::process(){
  std::string base_name = basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;
  ProteoformPtrVec proteo_ptrs 
      = readFastaToProteoform(db_file_name_, ResidueFactory::getBaseResiduePtrVec());
  PrsmPtrVec prsm_ptrs = readPrsm(input_file_name,proteo_ptrs);
  std::string output_file_name = base_name+"."+output_file_ext_;
  int max_id = prsm_ptrs[prsm_ptrs.size()-1]->getSpectrumId();
  PrsmWriter writer(output_file_name);
  for(int i=0; i<= max_id;i++){
    PrsmPtrVec selected_prsm_ptrs;
    for(size_t j=0;j<prsm_ptrs.size();j++){
      if(prsm_ptrs[j]->getSpectrumId()==i){
        selected_prsm_ptrs.push_back(prsm_ptrs[j]);
      }
    }
    PrsmPtrVec result_ptrs = getTopPrsms(selected_prsm_ptrs, n_top_);
    writer.writeVector(result_ptrs);
  }
  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  writer.close();
}
*/

void PrsmSelector::process(){
  std::string base_name = basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;
  PrsmReaderPtr reader_ptr(new PrsmReader(input_file_name));
  PrsmStrPtr prsm_str_ptr = reader_ptr->readOnePrsmStr();

  PrsmWriter writer(base_name +"."+output_file_ext_);
  
  // combine
  int spec_id = 0;
  while (prsm_str_ptr != nullptr) {
    PrsmStrPtrVec cur_str_ptrs;
    while (prsm_str_ptr != nullptr && prsm_str_ptr->getSpectrumId() == spec_id) {
      cur_str_ptrs.push_back(prsm_str_ptr);
      prsm_str_ptr = reader_ptr->readOnePrsmStr();
    }
    if (cur_str_ptrs.size() > 0) {
      std::sort(cur_str_ptrs.begin(),cur_str_ptrs.end(),prsmStrMatchEValueUp);
      for (int i = 0; i < top_num_; i++) {
        if (i >= (int)cur_str_ptrs.size()) {
          break;
        }
        writer.write(cur_str_ptrs[i]);
      }
    }
    spec_id++;
  }
  
  reader_ptr->close();
  writer.close();
}

} /* namespace prot */
