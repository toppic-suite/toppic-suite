#include "base/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_fdr.hpp"

namespace prot {
PrsmFdr::PrsmFdr(const std::string &db_file_name,
                 const std::string &spec_file_name,
                 const std::string &input_file_ext,
                 const std::string &output_file_ext): 
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext) {
    }

void PrsmFdr::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;

  PrsmStrPtrVec prsm_str_ptrs = PrsmReader::readAllPrsmStrs(input_file_name);

  PrsmStrPtrVec target_ptrs;
  PrsmStrPtrVec decoy_ptrs;
  for(size_t i=0; i< prsm_str_ptrs.size(); i++){
    if (prsm_str_ptrs[i]->getEValue() == 0.0) {
      LOG_ERROR("prot::PRSMFdr zero E value is reported");
    }
    else {
      std::string seq_name  = prsm_str_ptrs[i]->getSeqName();
      //LOG_DEBUG("seq name " << seq_name);
      if(seq_name.find("DECOY_")==0){
        decoy_ptrs.push_back(prsm_str_ptrs[i]);
      }
      else{
        target_ptrs.push_back(prsm_str_ptrs[i]);
      }
    }
  }
  compute(target_ptrs,decoy_ptrs);
  std::string output_file_name = base_name+"."+output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  std::sort(target_ptrs.begin(),target_ptrs.end(),PrsmStr::cmpSpectrumIdInc);
  writer.writeVector(target_ptrs);
  writer.close();
}

void PrsmFdr::compute(PrsmStrPtrVec &target_ptrs,PrsmStrPtrVec &decoy_ptrs){
  std::sort(target_ptrs.begin(),target_ptrs.end(),PrsmStr::cmpEValueInc);
  std::sort(decoy_ptrs.begin(),decoy_ptrs.end(),PrsmStr::cmpEValueInc);
  for(size_t i=0; i<target_ptrs.size(); i++){
    int n_target=i+1;
    double target_evalue = target_ptrs[i]->getEValue();
    int n_decoy = 0;
    for(size_t j=0; j<decoy_ptrs.size(); j++){
      if(decoy_ptrs[j]->getEValue() <= target_evalue){
        n_decoy++;
      }
      else{
        break;
      }
    }
    double fdr = (double)n_decoy/(double)n_target;
    if(fdr>1){
      fdr=1.0;
    }
    target_ptrs[i]->setFdr(fdr);
  }
}

} /* namespace prot */
