#include "base/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_top_selector.hpp"

namespace prot {

PrsmTopSelector::PrsmTopSelector(const std::string &db_file_name,
                                 const std::string &spec_file_name,
                                 const std::string &in_file_ext,
                                 const std::string &out_file_ext, 
                                 int n_top){
  spec_file_name_ = spec_file_name;
  db_file_name_ = db_file_name;
  input_file_ext_ = in_file_ext;
  output_file_ext_ = out_file_ext;
  n_top_ = n_top;
}

bool containsSameDbSeq(const PrsmStrPtrVec prsm_ptrs, PrsmStrPtr target_prsm_ptr) {
  for(size_t i=0; i< prsm_ptrs.size();i++){
    if(prsm_ptrs[i]->getDbSeqId() == target_prsm_ptr->getDbSeqId()){
      return true;
    }
  }
  return false;
}

PrsmStrPtrVec getTopPrsms(PrsmStrPtrVec &prsm_str_ptrs, int n_top){
  std::sort(prsm_str_ptrs.begin(),prsm_str_ptrs.end(),prsmStrEValueUp);
  int size = prsm_str_ptrs.size();
  int max = size > n_top? n_top:size;
  PrsmStrPtrVec result_ptrs;
  for(int i=0;i<max;i++){
    if(!containsSameDbSeq(result_ptrs, prsm_str_ptrs[i])){
      result_ptrs.push_back(prsm_str_ptrs[i]);
    }
  }
  return result_ptrs;
}

void PrsmTopSelector::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;
  PrsmReader reader(input_file_name);
  PrsmStrPtr prsm_str_ptr = reader.readOnePrsmStr();

  PrsmWriter writer(base_name +"."+output_file_ext_);
  
  int spec_id = 0;
  while (prsm_str_ptr != nullptr) {
    PrsmStrPtrVec cur_str_ptrs;
    while (prsm_str_ptr != nullptr && prsm_str_ptr->getSpectrumId() == spec_id) {
      cur_str_ptrs.push_back(prsm_str_ptr);
      prsm_str_ptr = reader.readOnePrsmStr();
    }
    PrsmStrPtrVec result_ptrs = getTopPrsms(cur_str_ptrs, n_top_);
    for (size_t i = 0; i < result_ptrs.size(); i++) {
      writer.write(result_ptrs[i]);
    }

    spec_id++;
  }
  
  reader.close();
  writer.close();
}

} /* namespace prot */
