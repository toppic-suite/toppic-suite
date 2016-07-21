#include "base/logger.hpp"

#include "spec/msalign_reader.hpp"
#include "feature/feature_detect_mng.hpp"
#include "feature/feature_detect.hpp"

namespace prot {

void readHeaders(std::string &file_name, MsHeaderPtrVec &header_ptr_vec) {
  MsAlignReader sp_reader(file_name, 1, nullptr);
  
  DeconvMsPtr ms_ptr;
  //LOG_DEBUG("Start search");
  while((ms_ptr = sp_reader.getNextMs())!= nullptr){
    header_ptr_vec.push_back(ms_ptr->getMsHeaderPtr());
    //std::cout << std::flush <<  "reading spectrum " << header_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  std::cout << std::endl;
}

void outputHeaders(MsHeaderPtrVec &header_ptr_vec) {
  for (size_t i = 0; i < header_ptr_vec.size(); i++) {
    MsHeaderPtr ptr = header_ptr_vec[i];
    std::cout << ptr->getId() << "\t" << ptr->getFirstScanNum() << "\t";
    std::cout <<  ptr->getRetentionTime() << "\t" << ptr->getPrecMonoMass() << "\t";
    std::cout << ptr->getPrecMonoMz() << "\t" << ptr->getPrecCharge() << "\t" << ptr->getPrecInte() << std::endl; 
  }

}

bool isConsistent(MsHeaderPtr &a, MsHeaderPtr &b, FeatureDetectMngPtr mng_ptr) {
  double mass_diff = std::abs(a->getPrecMonoMass() - b->getPrecMonoMass());
  double error_tole = a->getPrecMonoMass() * mng_ptr->ppo_;
  if (mass_diff <= error_tole) {
    return true;
  }
  return false;
}

void groupHeaders(MsHeaderPtrVec &header_ptr_vec, FeatureDetectMngPtr mng_ptr) {
  int total_num = header_ptr_vec.size();
  MsHeaderPtrVec remain_ptrs = header_ptr_vec; 
  MsHeaderPtrVec sorted_ptrs = remain_ptrs; 
  MsHeaderPtr2D results;
  while (sorted_ptrs.size() > 0) {
    std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), MsHeader::cmpPrecInteDec);
    MsHeaderPtr best_ptr = sorted_ptrs[0];
    MsHeaderPtrVec cur_group;
    cur_group.push_back(best_ptr);
    int best_id = best_ptr->getId();
    remain_ptrs[best_id] = nullptr;
    int intv_bgn = best_id - mng_ptr->intv_width_;
    if (intv_bgn < 0) {
      intv_bgn = 0;
    }
    int intv_end = best_id + mng_ptr->intv_width_;
    if (intv_end >= total_num) {
      intv_end = total_num -1;
    }
    for (int i = intv_bgn; i <= intv_end; i++) {
      if (remain_ptrs[i] != nullptr) {
        if (isConsistent(best_ptr, remain_ptrs[i], mng_ptr)) {
          cur_group.push_back(remain_ptrs[i]);
          remain_ptrs[i] = nullptr;
        }
      }
    }
    results.push_back(cur_group);
    sorted_ptrs.clear();
    for (size_t i = 0; i < remain_ptrs.size(); i++) {
      if (remain_ptrs[i] != nullptr) {
        sorted_ptrs.push_back(remain_ptrs[i]);
      }
    }
    LOG_DEBUG("sorted ptr size " << sorted_ptrs.size());
  }
  for (size_t i = 0; i < results.size(); i++) {
    std::cout << "Group " << i << " number " << results[i].size() << std::endl;
    for (size_t j = 0; j < results[i].size(); j++) {
      MsHeaderPtr ptr = results[i][j];
      std::cout << "\t" << ptr->getId() << "\t" << ptr->getFirstScanNum() << "\t";
      std::cout <<  ptr->getRetentionTime() << "\t" << ptr->getPrecMonoMass() << "\t";
      std::cout << ptr->getPrecMonoMz() << "\t" << ptr->getPrecCharge() << "\t" << ptr->getPrecInte() << std::endl; 
    }
  }

}

void FeatureDetect::process(std::string &input_file_name){
  FeatureDetectMngPtr mng_ptr(new FeatureDetectMng());
  MsHeaderPtrVec header_ptr_vec;
  readHeaders(input_file_name, header_ptr_vec);
  groupHeaders(header_ptr_vec, mng_ptr);
  //outputHeaders(header_ptr_vec);
}

//std::string output_file_name = FileUtil::basename(sp_file_name)+".feature";
}
