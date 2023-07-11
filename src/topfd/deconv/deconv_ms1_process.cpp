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
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "ms/spec/msalign_writer.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"

#include "topfd/deconv/deconv_one_sp.hpp"
#include "topfd/deconv/deconv_ms1_process.hpp"

/*
#include "common/util/version.hpp"
#include "common/util/time_util.hpp"

#include "ms/spec/baseline_util.hpp"
#include "ms/spec/msalign_thread_merge.hpp"

#include "ms/env/env_base.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_writer.hpp"

#include "ms/mzml/mzml_ms_json_writer.hpp"

#include "topfd/common/topfd_para.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp" 
*/

namespace toppic {

DeconvMs1Process::DeconvMs1Process(TopfdParaPtr topfd_para_ptr) {
  topfd_para_ptr_ = topfd_para_ptr;
}

std::string DeconvMs1Process::updateMsg(MsHeaderPtr header_ptr, int scan_cnt, 
                                        int total_scan_num) {
  std::string percentage = str_util::toString(scan_cnt * 100 / total_scan_num);
  std::string msg = "Processing spectrum scan " 
    + std::to_string(header_ptr->getFirstScanNum()) + " ...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

void DeconvMs1Process::prepareFileFolder() {
  if (topfd_para_ptr_->isGeneHtmlFolder()) {
    //json file names
    std::string html_dir = topfd_para_ptr_->getHtmlDir();
    if (!file_util::exists(html_dir)) {
      file_util::createFolder(html_dir);
    }
    std::string ms1_json_dir = topfd_para_ptr_->getMs1JsonDir();
    if (!file_util::exists(ms1_json_dir)) {
      file_util::createFolder(ms1_json_dir);
    }
  }
}

std::function<void()> geneTask(MzmlMsGroupPtr ms_group_ptr, 
                               TopfdParaPtr topfd_para_ptr, 
                               MsAlignWriterPtrVec ms1_writer_ptr_vec, 
                               SimpleThreadPoolPtr pool_ptr) { 
  return [ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr]() {

    //deconvMsOne(ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr);

  };
}


void DeconvMs1Process::process() {
  MzmlMsGroupReaderPtr reader_ptr = 
    std::make_shared<MzmlMsGroupReader>(mzml_file_name_, 
                                        topfd_para_ptr_->getPrecWindow(),  
                                        topfd_para_ptr_->getActivation(),
                                        topfd_para_ptr_->getFracId(),
                                        topfd_para_ptr_->isFaims(), 
                                        topfd_para_ptr_->getFaimsVoltage(), 
                                        topfd_para_ptr_->isMissingLevelOne());

  MzmlMsGroupPtr ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  if (ms_group_ptr == nullptr) {
    LOG_ERROR("No spectrum to read in mzML file!");
    return;
  }
  prepareFileFolder();
  // init thread pool
  int thread_num = topfd_para_ptr_->getThreadNum();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);  
  // init msalign writer vector for multiple threads
  std::string ms1_msalign_name = output_base_name_ + "_ms1.msalign";
  MsAlignWriterPtrVec ms1_writer_ptr_vec;
  for (int i = 0; i < thread_num; i++) { 
    MsAlignWriterPtr ms1_ptr 
        = std::make_shared<MsAlignWriter>(ms1_msalign_name + "_" + str_util::toString(i));
    ms1_writer_ptr_vec.push_back(ms1_ptr);
  }
  // counter for processed spectra
  int spec_cnt = 0;
  // total spectrum number
  int total_spec_num = topfd_para_ptr_->getMs1ScanNum(); 
  while (ms_group_ptr != nullptr) {
    while(pool_ptr->getQueueSize() >= thread_num * 2){
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneTask(ms_group_ptr, topfd_para_ptr_, ms1_writer_ptr_vec, pool_ptr)); 
    spec_cnt++; 
    std::string msg = updateMsg(ms_group_ptr->getMsOnePtr()->getMsHeaderPtr(), 
                                spec_cnt, total_spec_num);
    std::cout << "\r" << msg << std::flush;
    ms_group_ptr = reader_ptr->getNextMsGroupPtr();    
  }
  pool_ptr->ShutDown();
  /*
  // Merge files
  std::string ms1_file_name = "ms1.msalign";
  MsalignThreadMergePtr ms1_merge_ptr
            = std::make_shared<MsalignThreadMerge>(spec_file_name_, ms1_file_name, 
                                                  thread_num_, ms1_file_name, 
                                                  topfd_para_ptr_-> getParaStr("#", "\t"));
  ms1_merge_ptr->process();
  // Remove temporary files to add
  //

  */

  std::cout << std::endl;
}


/*

void deconvMsOne(MzmlMsGroupPtr ms_group_ptr, 
                 TopfdParaPtr topfd_para_ptr,  
                 MsAlignWriterPtrVec ms1_writer_ptr_vec, 
                 SimpleThreadPoolPtr pool_ptr) { 
  // 1. Store peak intensity 
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  std::vector<double> intensities;
  for (size_t i = 0; i < peak_list.size(); i++) {
    intensities.push_back(peak_list[i]->getIntensity());
  }
  // 2. Deconv envelopes in precursor windows and remove them
  MatchEnvPtrVec prec_envs = RawMsGroupFaimeReader::obtainPrecEnvs(ms_group_ptr, 
                                                                   topfd_para_ptr->getMaxMass(),
                                                                   topfd_para_ptr->getMaxCharge()); 
  //remove precursor peaks
  for (size_t i = 0; i < prec_envs.size(); i++) {
    ExpEnvPtr env_ptr = prec_envs[i]->getRealEnvPtr();
    for (int p = 0; p < env_ptr->getPeakNum(); p++) {
      if (env_ptr->isExist(p)) {
        int idx = env_ptr->getPeakIdx(p);
        peak_list[idx]->setIntensity(0);
      }
    }
  }
  // 3. Deconv the whole spectrum 
  MatchEnvPtrVec deconv_envs;
  if (peak_list.size() > 0) {
    LOG_DEBUG("set data....");
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    DeconvSingleSpPtr deconv_ptr = std::make_shared<DeconvSingleSp>(topfd_para_ptr_, peak_list, ms_level);
    deconv_ptr->run();
    deconv_envs = deconv_ptr->getResult();
  }
  // 4. Filter only deconvolution envelopes, precursor envelopes are not
  // filtered.
  
  // 5. Merge precursor envelopes and deconvolution envelopes
  MatchEnvPtrVec result_envs;
  result_envs.insert(result_envs.end(), prec_envs.begin(), prec_envs.end());
  result_envs.insert(result_envs.end(), deconv_envs.begin(), deconv_envs.end());
  LOG_DEBUG("result num " << prec_envs.size());
  
  // 6. Write to msalign file
  DeconvMsPtr deconv_ms_ptr = match_env_util::getDeconvMsPtr(header_ptr,
                                                             result_envs,
                                                             topfd_para_ptr->isUseEnvCnn());

  boost::thread::id thread_id = boost::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);
  ms1_writer_ptr_vec[writer_id]->writeMs(deconv_ms_ptr);
  
  //7. write json file 
  if (gene_html_dir) {
    // add back precursor peaks
    for (size_t i = 0; i < peak_list.size(); i++) {
      peak_list[i]->setIntensity(intensities[i]);
    }
    std::string json_file_name = ms1_json_dir 
        + file_util::getFileSeparator() 
        + "spectrum" + std::to_string(header_ptr->getSpecId()) + ".js";
    mzml_ms_json_writer::write(json_file_name, ms_ptr, prec_envs);    
  }
}
*/

}; // namespace toppic
