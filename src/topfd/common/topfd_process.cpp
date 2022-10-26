//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/time_util.hpp"
#include "common/base/base_data.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/spec/deconv_json_merge.hpp"
#include "ms/feature/feature_merge.hpp"
#include "ms/env/env_base.hpp"
#include "topfd/deconv/deconv_process.hpp"
#include "topfd/feature_detect/feature_detect.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

namespace topfd_process {

void processOneFile(TopfdParaPtr para_ptr,  
                    const std::string &spec_file_name, 
                    int frac_id) {
  try {
    int thread_number = para_ptr->getThreadNum();

    std::cout << "Processing " << spec_file_name << " started." << std::endl;
    std::cout << "Deconvolution started." << std::endl;

    DeconvProcess processor(para_ptr, spec_file_name, frac_id, thread_number);
    
    processor.process();
    std::cout << "Deconvolution finished." << std::endl;

    std::cout << "Deleting temporary files - started." << std::endl;
    
    std::string file_num = "";

    for (int i = 0; i < processor.msalign_num_; i++) {
      if (processor.isFaims_) {file_num = str_util::toString(i) + "_";} // if FAIME Data

      std::string base_name_ms1 = file_util::basename(spec_file_name) + "_" + file_num + "ms1.msalign";
      std::string base_name_ms2 = file_util::basename(spec_file_name) + "_" + file_num + "ms2.msalign";

      std::string fa_base_ms1 = file_util::absoluteName(base_name_ms1);
      std::string fa_base_ms2 = file_util::absoluteName(base_name_ms2);

      std::replace(fa_base_ms1.begin(), fa_base_ms1.end(), '\\', '/');
      std::replace(fa_base_ms2.begin(), fa_base_ms2.end(), '\\', '/');

      file_util::cleanPrefix(base_name_ms1 + "_", fa_base_ms1 + "_");
      file_util::cleanPrefix(base_name_ms2 + "_", fa_base_ms2 + "_");
    }
    std::cout << "Deleting temporary files - finished." << std::endl; 

    std::cout << "Feature detection started." << std::endl;
    feature_detect::process(frac_id, 
                            spec_file_name,
                            para_ptr->isMissingLevelOne(), 
                            para_ptr->getResourceDir(),             
			                      para_ptr->getActivation(), 
                            processor.isFaims_,
                            processor.voltage_vec_);
    std::cout << "Feature detection finished." << std::endl;
    
    std::cout << "Processing " << spec_file_name << " finished." << std::endl;
  } catch (const char* e) {
    LOG_ERROR("[Exception] " << e);
    exit(EXIT_FAILURE);
  } 
}

void moveFiles(std::string &spec_file_name, bool move_mzrt) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string file_dir =  base_name + "_file";
  std::string file_name = base_name + "_ms1.msalign";

  //create folder only if ms1 msalign and frac mzrt csv exist  
  //== when ms1 spectra was used
  //non-FAIME data
  if (file_util::exists(file_name)) {//if there is ms1 msalign
    if (!file_util::exists(file_dir)) {
      file_util::createFolder(file_dir);
    }
    file_util::moveFile(file_name, file_dir);
    //file_name = base_name + "_feature.xml";
    //file_util::moveFile(file_name, file_dir);
    if (move_mzrt) {
      file_name = base_name + "_frac.mzrt.csv";
      if (file_util::exists(file_name)) {
        file_util::moveFile(file_name, file_dir);
      }
    }
  }
  //FAIME data
  int file_num = 0; 

  file_dir =  base_name + "_" + str_util::toString(file_num) + "_file";
  file_name = base_name + "_" + str_util::toString(file_num) + "_ms1.msalign";
  std::string frac_file_name = base_name + "_" + str_util::toString(file_num) + "_frac.mzrt.csv";

  while (file_util::exists(file_name)) {
    if (!file_util::exists(file_dir)) {
      file_util::createFolder(file_dir);
    }
    file_util::moveFile(file_name, file_dir);

    if (move_mzrt) {
      file_util::moveFile(frac_file_name, file_dir);
    }
    file_num++;
    file_name = base_name + "_" + str_util::toString(file_num) + "_ms1.msalign";
    file_dir =  base_name + "_" + str_util::toString(file_num) + "_file";
    frac_file_name = base_name + "_" + str_util::toString(file_num) + "_frac.mzrt.csv";
  }  
    //file_name = base_name + "_feature.xml";
    //file_util::moveFile(file_name, file_dir);
  /*
  if (move_sample_feature) {
    file_name = base_name + "_ms1.feature";
    file_util::moveFile(file_name, file_dir);
  }
  */
}

bool isValidFile(std::string &file_name) {
  if (str_util::endsWith(file_name, "mzML")
      || str_util::endsWith(file_name, "mzXML")
      || str_util::endsWith(file_name, "mzml")
      || str_util::endsWith(file_name, "mzxml")) {
    return true;
  }
  else {
    return false;
  }
}

/*
void mergeFiles(TopfdParaPtr para_ptr, 
                std::vector<std::string> &spec_file_lst) {
  try{
    std::string para_str = para_ptr->getParaStr("#");
    std::string merged_file_name = para_ptr->merged_file_name_;
    std::cout << "Merging files started." << std::endl;
    MsAlignFracMergePtr msalign_merger 
        = std::make_shared<MsAlignFracMerge>(spec_file_lst, merged_file_name);
    msalign_merger->process(para_str);
    msalign_merger = nullptr;
    DeconvJsonMergePtr json_merger 
        = std::make_shared<DeconvJsonMerge>(spec_file_lst, merged_file_name);
    json_merger->process();
    json_merger = nullptr;
    feature_merge::process(spec_file_lst, merged_file_name, para_str);
    std::cout << "Merging files finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
}
*/

void getMsScanCount(std::string spectrum_file_name, std::vector<int> &scan_cnt_vec) {
  typedef std::shared_ptr<pwiz::msdata::MSDataFile> MSDataFilePtr;
  pwiz::msdata::DefaultReaderList readers_;
  MSDataFilePtr msd_ptr_ = std::make_shared<pwiz::msdata::MSDataFile>(spectrum_file_name, &readers_);
  pwiz::msdata::SpectrumListPtr spec_list_ptr_ =  msd_ptr_->run.spectrumListPtr;
  
  int sp_size = spec_list_ptr_->size();
  int ms1_cnt = 0;
  int ms2_cnt = 0;

  for(int i = 0; i < sp_size; i++){
      int scan_level = spec_list_ptr_->spectrum(i)->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>(); // check scanLevel
      if (scan_level == 2) {
        ms2_cnt++;
      }
      else if (scan_level == 1) {
        ms1_cnt++;
      }
  }
  scan_cnt_vec.push_back(ms1_cnt);
  scan_cnt_vec.push_back(ms2_cnt);
}

int process(TopfdParaPtr para_ptr,  std::vector<std::string> spec_file_lst) {
  base_data::init();

  EnvBase::initBase(para_ptr->getResourceDir());
  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    //get ms1 and ms2 scan number
    std::vector<int> scan_cnt_vec;
    getMsScanCount(spec_file_lst[k], scan_cnt_vec);

    para_ptr->setMs1ScanNumber(scan_cnt_vec[0]);
    para_ptr->setMs2ScanNumber(scan_cnt_vec[1]);

    std::string print_str = para_ptr->getParaStr("", " ");//print parameter for each file
    std::cout << print_str;

    if (isValidFile(spec_file_lst[k])) {
      processOneFile(para_ptr, spec_file_lst[k], k);
      bool move_mzrt = true;
      moveFiles(spec_file_lst[k], move_mzrt); 
      std::cout << "Timestamp: " << time_util::getTimeStr() << std::endl;
    }
    else {
      std::cout << spec_file_lst[k] << " is not a valid mass spectral file!" << std::endl; 
    }
  }

  base_data::release();
  std::cout << "TopFD finished." << std::endl << std::flush;

  return 0;
}

} // namespace topfd_process 

}  // namespace toppic
