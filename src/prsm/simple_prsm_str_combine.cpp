#include "base/file_util.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace prot {

SimplePrsmStrCombine::SimplePrsmStrCombine(const std::string &spec_file_name, 
                                           const std::vector<std::string> &in_file_exts,
                                           const std::string &out_file_ext, 
                                           int top_num) {
  input_file_exts_ = in_file_exts;
  output_file_ext_ = out_file_ext;
  spec_file_name_ = spec_file_name;
  top_num_ = top_num;
}

SimplePrsmStrCombine::SimplePrsmStrCombine(const std::string &spec_file_name, 
                                           const std::string &in_file_ext,
                                           int in_num,
                                           const std::string &out_file_ext, 
                                           int top_num) {
  output_file_ext_ = out_file_ext;
  spec_file_name_ = spec_file_name;
  top_num_ = top_num;
  for (int i = 0; i < in_num; i ++) {
    std::string ext = in_file_ext + "_" + std::to_string(i);
    input_file_exts_.push_back(ext);
  }
}

void SimplePrsmStrCombine::process() {
  size_t input_num = input_file_exts_.size();
  std::string base_name = basename(spec_file_name_); 
  // open files
  SimplePrsmReaderPtrVec reader_ptrs;
  SimplePrsmStrPtrVec prsm_str_ptrs;
  for (size_t i = 0; i < input_num; i++) {
    std::string input_file_name = base_name + "." + input_file_exts_[i]; 
    SimplePrsmReaderPtr reader_ptr(new SimplePrsmReader(input_file_name));
    LOG_DEBUG("input file name " << input_file_name);
    SimplePrsmStrPtr str_ptr = reader_ptr->readOnePrsmStr();
    reader_ptrs.push_back(reader_ptr);
    prsm_str_ptrs.push_back(str_ptr);
  }
  SimplePrsmWriter writer(base_name +"."+output_file_ext_);
  
  // combine
  int spec_id = 0;
  bool finish = false;
  while (!finish) {
    //LOG_DEBUG("spec id " << spec_id);
    finish = true;
    SimplePrsmStrPtrVec cur_str_ptrs;
    for (size_t i = 0; i < input_num; i++) {
      if (prsm_str_ptrs[i] != nullptr) {
        finish = false;
        while (prsm_str_ptrs[i]!= nullptr && 
               prsm_str_ptrs[i]->getSpectrumId() == spec_id) {
          cur_str_ptrs.push_back(prsm_str_ptrs[i]);
          prsm_str_ptrs[i] = reader_ptrs[i]->readOnePrsmStr();
        }
      }
    }
    if (cur_str_ptrs.size() > 0) {
      std::sort(cur_str_ptrs.begin(),cur_str_ptrs.end(),simplePrsmStrScoreDown);
      for (int i = 0; i < top_num_; i++) {
        if (i >= (int)cur_str_ptrs.size()) {
          break;
        }
        writer.write(cur_str_ptrs[i]);
      }
    }
    spec_id++;
  }
  
  // close files
  for (size_t i = 0; i < input_num; i++) {
    reader_ptrs[i]->close();
  }

  writer.close();
}

} /* namespace prot */
