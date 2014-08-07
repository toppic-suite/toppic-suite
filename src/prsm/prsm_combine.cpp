#include "base/file_util.hpp"
#include "prsm/prsm_combine.hpp"

namespace prot {

PrsmCombine::PrsmCombine(const std::string &db_file_name, 
                         const std::string &spec_file_name, 
                         const std::vector<std::string> &in_file_exts,
                         const std::string &out_file_ext) {
  input_file_exts_ = in_file_exts;
  output_file_ext_ = out_file_ext;
  spec_file_name_ = spec_file_name;
  db_file_name_ = db_file_name;
}

void PrsmCombine::process() {
  ProteoformPtrVec proteo_ptrs 
      = readFastaToProteoform(db_file_name_, ResidueFactory::getBaseResiduePtrVec());
  PrsmPtrVec prsm_ptrs;
  for(size_t i=0;i<input_file_exts_.size();i++){
    PrsmPtrVec temp_ptrs 
        = readPrsm(basename(spec_file_name_)+"."+input_file_exts_[i],proteo_ptrs);
    for(size_t j=0;j<temp_ptrs.size();j++){
      prsm_ptrs.push_back(temp_ptrs[j]);
    }
  }
  std::sort(prsm_ptrs.begin(),prsm_ptrs.end(),prsmSpectrumIdUpMatchFragDown);
  PrsmWriter all_writer(basename(spec_file_name_)+"."+output_file_ext_);
  all_writer.writeVector(prsm_ptrs);

  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  all_writer.close();
}

} /* namespace prot */
