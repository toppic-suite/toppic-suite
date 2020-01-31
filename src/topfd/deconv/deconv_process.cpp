//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "ms/spec/baseline_util.hpp"
#include "common/util/version.hpp"
#include "common/util/time_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_writer.hpp"

#include "topfd/common/topfd_para.hpp"
#include "topfd/msreader/raw_ms_writer.hpp"
#include "topfd/deconv/deconv_process.hpp"

#include <mutex>
#include <vector>

namespace toppic {

std::mutex count_lock;
std::mutex ms1_map_lock;
std::mutex ms2_map_lock;

DeconvProcess::DeconvProcess(TopfdParaPtr topfd_para_ptr, 
                             const std::string &spec_file_name, 
                             int frac_id, int t, DeconvProcess *deconv) {
  topfd_para_ptr_ = topfd_para_ptr;
  env_para_ptr_ = std::make_shared<EnvPara>(topfd_para_ptr);
  dp_para_ptr_ = std::make_shared<DpPara>();
  deconv_instance_ptr = deconv;

  spec_file_name_ = spec_file_name;
  frac_id_ = frac_id; 
  thread_num = t;
  prepareFileFolder();

}

void DeconvProcess::writeMsalign(MsAlignWriterPtr ms_wtr_ptr, std::map<int, DeconvMsPtr> ms_ptr_map){
  for (int i = 0; i < ms_ptr_map.size(); i++){
    ms_wtr_ptr->write(ms_ptr_map[i]);
  }
}

std::string DeconvProcess::updateMsg(MsHeaderPtr header_ptr, int scan, 
                                     int total_scan_num) {
  std::string percentage = str_util::toString(scan * 100 / total_scan_num);
  std::string msg = "Processing spectrum scan " 
      + std::to_string(header_ptr->getFirstScanNum()) + " ...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

void DeconvProcess::prepareFileFolder() {
  base_name_ = file_util::basename(spec_file_name_);

  //envelope file names
  ms1_env_name_ = base_name_ + "_ms1.env"; 
  ms2_env_name_ = base_name_ + "_ms2.env"; 
  if (file_util::exists(ms1_env_name_)) {
    file_util::delFile(ms1_env_name_); 
  }
  if (file_util::exists(ms2_env_name_)) {
    file_util::delFile(ms2_env_name_); 
  }

  //json file names
  html_dir_ =  base_name_ + "_html";
  ms1_json_dir_ = html_dir_ 
      + file_util::getFileSeparator() + "topfd" 
      + file_util::getFileSeparator() + "ms1_json";
  ms2_json_dir_ = html_dir_ 
      + file_util::getFileSeparator() + "topfd" 
      + file_util::getFileSeparator() + "ms2_json";
  if (file_util::exists(html_dir_)) {
    file_util::delDir(html_dir_);
  }
  file_util::createFolder(html_dir_);
  file_util::createFolder(ms1_json_dir_);
  file_util::createFolder(ms2_json_dir_);
}

void DeconvProcess::process() {
  // writer 
  std::string ms1_msalign_name, ms2_msalign_name;
  ms1_msalign_name = base_name_ + "_ms1.msalign";
  ms2_msalign_name = base_name_ + "_ms2.msalign";

  MsAlignWriterPtr ms1_writer_ptr = std::make_shared<MsAlignWriter>(ms1_msalign_name);
  MsAlignWriterPtr ms2_writer_ptr = std::make_shared<MsAlignWriter>(ms2_msalign_name);

  std::string para_str = topfd_para_ptr_->getParaStr("#");

  ms1_writer_ptr->writePara(para_str);
  ms2_writer_ptr->writePara(para_str);

  // reader
  RawMsGroupReaderPtr reader_ptr = std::make_shared<RawMsGroupReader>(spec_file_name_, 
                                                                      topfd_para_ptr_->missing_level_one_,
                                                                      frac_id_);
  RawMsGroupReaderPtr reader_pl_ptr = std::make_shared<RawMsGroupReader>(spec_file_name_, 
                                                                      topfd_para_ptr_->missing_level_one_,
                                                                      frac_id_);

  // processor
  DeconvOneSpPtr deconv_ptr = std::make_shared<DeconvOneSp>(env_para_ptr_, dp_para_ptr_);

  if (topfd_para_ptr_->missing_level_one_) {
    processSpMissingLevelOne(deconv_ptr, reader_ptr, ms2_writer_ptr);
  }
  else {
    processSp(deconv_ptr, reader_ptr, reader_pl_ptr, ms1_writer_ptr, ms2_writer_ptr);
  }
}
void DeconvProcess::deconvMissingMsOne(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, std::map<int, DeconvMsPtr> *ms1_missing_map_ptr){
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();

  header_ptr->setPrecCharge(EnvPara::getDefaultMaxCharge());
  double prec_mz = EnvPara::getDefaultMaxMass()/EnvPara::getDefaultMaxCharge();
  header_ptr->setPrecMonoMz(prec_mz);
  header_ptr->setPrecSpMz(prec_mz);
  MatchEnvPtrVec result_envs; 

  if (peak_list.size() > 0) {
    deconv_ptr->setData(peak_list, EnvPara::getDefaultMaxMass(),
                            EnvPara::getDefaultMaxCharge());
    deconv_ptr->run();
    result_envs = deconv_ptr->getResult();
  }

  DeconvMsPtr deconv_ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);

  int index = header_ptr->getId();

  ms2_map_lock.lock();
  (*ms1_missing_map_ptr).insert ( std::pair<int, DeconvMsPtr>(index,deconv_ms_ptr) );
  ms2_map_lock.unlock();

  if (topfd_para_ptr_->output_match_env_) {
    match_env_writer::write(ms2_env_name_, header_ptr, result_envs);
  }

  if (topfd_para_ptr_->output_json_files_) {
    std::string json_file_name = ms2_json_dir_ 
        + file_util::getFileSeparator() 
        + "spectrum" + std::to_string(header_ptr->getId()) + ".js";
    raw_ms_writer::write(json_file_name, ms_ptr, result_envs);    
  }

}
std::function<void()> geneTaskMissingMs1(RawMsGroupPtr ms_group_ptr, DeconvOneSpPtr deconv_ptr,
                                DeconvProcess *deconv_instance_ptr, int **count2, int total_scan_num, std::map<int, DeconvMsPtr> *ms1_missing_map_ptr){
  return [ms_group_ptr, deconv_ptr, deconv_instance_ptr, count2, total_scan_num, ms1_missing_map_ptr]() {
  
  RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();                           
  
  for (size_t i = 0; i < ms_two_ptr_vec.size(); i++) {
      RawMsPtr ms_two_ptr = ms_two_ptr_vec[i];
      std::string msg = deconv_instance_ptr->updateMsg(ms_two_ptr->getMsHeaderPtr(), **count2 + 1, total_scan_num);
      std::stringstream msgStream;
      msgStream << std::flush << msg << "\r";
      
      deconv_instance_ptr->deconvMissingMsOne(ms_two_ptr, deconv_ptr, ms1_missing_map_ptr);

      count_lock.lock();

      std::cout << msgStream.str();
      **count2 = **count2 + 1;

      count_lock.unlock();
  }
  };
}

void DeconvProcess::processSpMissingLevelOne(DeconvOneSpPtr deconv_ptr, 
                                             RawMsGroupReaderPtr reader_ptr,
                                             MsAlignWriterPtr ms2_writer_ptr) {
  // reader_ptr
  int total_scan_num = reader_ptr->getInputSpNum();
  RawMsGroupPtr ms_group_ptr;
  ms_group_ptr = reader_ptr->getNextMsGroupPtr();

  int count2 = 0;
  int *count2_ptr;
  count2_ptr = &count2;

  std::map<int, DeconvMsPtr> ms1_missing_map; 
  std::map<int, DeconvMsPtr> *ms1_missing_map_ptr;
  ms1_missing_map_ptr = &ms1_missing_map;

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);  

  while (ms_group_ptr != nullptr) {
    while(pool_ptr->getQueueSize() >= thread_num * 2){
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }

    //to generate a new deconOneSp object each time
    EnvParaPtr env_ptr = deconv_ptr->getEnvPara();
    DpParaPtr dp_ptr = deconv_ptr->getDpPara();

    std::shared_ptr<EnvPara> env_ptr_new = std::make_shared<EnvPara>(*env_ptr);
    std::shared_ptr<DpPara> dp_ptr_new = std::make_shared<DpPara>(*dp_ptr);

    DeconvOneSpPtr deconv_ptr_new = std::make_shared<DeconvOneSp>(env_ptr_new, dp_ptr_new);
    
    pool_ptr->Enqueue(geneTaskMissingMs1(ms_group_ptr, deconv_ptr_new, deconv_instance_ptr, &count2_ptr, total_scan_num, ms1_missing_map_ptr));
    ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  }
  pool_ptr->ShutDown();

  writeMsalign(ms2_writer_ptr, ms1_missing_map);

  std::cout << std::endl;
}

void DeconvProcess::deconvMsOne(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, 
                                MatchEnvPtrVec prec_envs, MsAlignWriterPtr ms1_writer_ptr, std::map<int, DeconvMsPtr> *ms1_ptr_map_ptr) { 
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  LOG_DEBUG("peak list size " << peak_list.size());
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  LOG_DEBUG("ms level " << header_ptr->getMsLevel() );
  // int scan_num_ = header_ptr->getFirstScanNum();
  //remove precursor peaks
  std::vector<double> intensities;
  for (size_t i = 0; i < peak_list.size(); i++) {
    intensities.push_back(peak_list[i]->getIntensity());
  }
  for (size_t i = 0; i < prec_envs.size(); i++) {
    RealEnvPtr env_ptr = prec_envs[i]->getRealEnvPtr();
    for (int p = 0; p < env_ptr->getPeakNum(); p++) {
      if (env_ptr->isExist(p)) {
        int idx = env_ptr->getPeakIdx(p);
        peak_list[idx]->setIntensity(0);
      }
    }
  }
  MatchEnvPtrVec result_envs;

  if (peak_list.size() > 0) {
    LOG_DEBUG("set data....");
    deconv_ptr->setMsLevel(header_ptr->getMsLevel());
    deconv_ptr->setData(peak_list);
    deconv_ptr->run();
    result_envs = deconv_ptr->getResult();
  }
  prec_envs.insert(prec_envs.end(), result_envs.begin(), result_envs.end());
  LOG_DEBUG("result num " << prec_envs.size());

  DeconvMsPtr deconv_ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, prec_envs);
  
  int index = header_ptr->getId();

  ms1_map_lock.lock();
  (*ms1_ptr_map_ptr).insert ( std::pair<int, DeconvMsPtr>(index,deconv_ms_ptr) );
  ms1_map_lock.unlock();

  if (topfd_para_ptr_->output_match_env_) {
    match_env_writer::write(ms1_env_name_, header_ptr, prec_envs);
  }

  // add back precursor peaks
  for (size_t i = 0; i < peak_list.size(); i++) {
    peak_list[i]->setIntensity(intensities[i]);
  }
  if (topfd_para_ptr_->output_json_files_) {
    std::string json_file_name = ms1_json_dir_ 
        + file_util::getFileSeparator() 
        + "spectrum" + std::to_string(header_ptr->getId()) + ".js";
    raw_ms_writer::write(json_file_name, ms_ptr, prec_envs);    
  }
}

void DeconvProcess::deconvMsTwo(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, 
                                MsAlignWriterPtr ms2_writer_ptr, std::map<int, DeconvMsPtr> *ms2_ptr_map_ptr) { 

  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  LOG_DEBUG("peak list size " << peak_list.size());
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  LOG_DEBUG("ms level " << header_ptr->getMsLevel() );
  // int scan_num_ = header_ptr->getFirstScanNum();
  double max_frag_mass = header_ptr->getPrecMonoMass();

  if (max_frag_mass == 0.0) {
    max_frag_mass = header_ptr->getPrecSpMass();
  }
    deconv_ptr->setMsLevel(header_ptr->getMsLevel());
    deconv_ptr->setData(peak_list, max_frag_mass, header_ptr->getPrecCharge());
    deconv_ptr->run();
    MatchEnvPtrVec result_envs = deconv_ptr->getResult();
    DeconvMsPtr deconv_ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);

    int index = header_ptr->getId();
    ms2_map_lock.lock();
    (*ms2_ptr_map_ptr).insert ( std::pair<int, DeconvMsPtr>(index,deconv_ms_ptr) );
     ms2_map_lock.unlock();
    //ms2_writer_ptr->write(deconv_ms_ptr);

    if (topfd_para_ptr_->output_match_env_) {
      match_env_writer::write(ms2_env_name_, header_ptr, result_envs);
    }
    if (topfd_para_ptr_->output_json_files_) {
      std::string json_file_name = ms2_json_dir_ 
          + file_util::getFileSeparator() 
          + "spectrum" + std::to_string(ms_ptr->getMsHeaderPtr()->getId()) + ".js";
      raw_ms_writer::write(json_file_name, ms_ptr, result_envs);    
    }

}

//DecovOne & Two
std::function<void()> geneTask(RawMsGroupPtr ms_group_ptr, DeconvOneSpPtr deconv_ptr, MatchEnvPtrVec prec_envs,
                                MsAlignWriterPtr ms_writer_ptr, MsAlignWriterPtr ms2_writer_ptr, 
                                DeconvProcess *deconv_instance_ptr, int **count1, int **count2, int total_scan_num, std::map<int, DeconvMsPtr> *ms1_ptr_map_ptr, std::map<int, DeconvMsPtr> *ms2_ptr_map_ptr){
  return [ms_group_ptr, deconv_ptr, prec_envs, ms_writer_ptr, ms2_writer_ptr, deconv_instance_ptr, count1, count2, total_scan_num, ms1_ptr_map_ptr, ms2_ptr_map_ptr]() {
   
   RawMsPtr ms_one_ptr = ms_group_ptr->getMsOnePtr();                            
   RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();

   std::string msg = deconv_instance_ptr->updateMsg(ms_one_ptr->getMsHeaderPtr(), **count1 + **count2 + 1, total_scan_num);
   
   std::stringstream msgStream;
   msgStream << std::flush << msg << "\r";

   deconv_instance_ptr->deconvMsOne(ms_one_ptr, deconv_ptr, prec_envs, ms_writer_ptr, ms1_ptr_map_ptr);
  
   count_lock.lock();
   
   std::cout << msgStream.str();
   **count1 = **count1 + 1;
   
   count_lock.unlock();

  for (size_t i = 0; i < ms_two_ptr_vec.size(); i++) {
      RawMsPtr ms_two_ptr = ms_two_ptr_vec[i];
      msg = deconv_instance_ptr->updateMsg(ms_two_ptr->getMsHeaderPtr(), **count1 + **count2 + 1, total_scan_num);
      
      msgStream << std::flush << msg << "\r";
      
      deconv_instance_ptr->deconvMsTwo(ms_two_ptr, deconv_ptr, ms2_writer_ptr, ms2_ptr_map_ptr);

      count_lock.lock();

      std::cout << msgStream.str();
      **count2 = **count2 + 1;

      count_lock.unlock();
  }
    
  };
}
void DeconvProcess::processSp(DeconvOneSpPtr deconv_ptr, 
                              RawMsGroupReaderPtr reader_ptr, RawMsGroupReaderPtr reader_pl_ptr, MsAlignWriterPtr ms1_writer_ptr, MsAlignWriterPtr ms2_writer_ptr) {

  int total_scan_num = reader_ptr->getInputSpNum();

  RawMsGroupPtr ms_group_ptr;

  ms_group_ptr = reader_ptr->getNextMsGroupPtr();

  int count1 = 0;
  int count2 = 0;

  int *count1_ptr;
  int *count2_ptr;

  count1_ptr = &count1;
  count2_ptr = &count2;

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);  

  //generate vector that contains msalign writing information

  std::map<int, DeconvMsPtr> ms1_ptr_map; 
  std::map<int, DeconvMsPtr> *ms1_ptr_map_ptr;
  ms1_ptr_map_ptr = &ms1_ptr_map;

  std::map<int, DeconvMsPtr> ms2_ptr_map; 
  std::map<int, DeconvMsPtr> *ms2_ptr_map_ptr;
  ms2_ptr_map_ptr = &ms2_ptr_map;

  //generate writer ptr

  std::vector<MsAlignWriterPtr> writer_vec_1;
  std::vector<MsAlignWriterPtr> writer_vec_2;

  for (int i = 0; i < 8; i++){//change if 8 files are too small
    MsAlignWriterPtr ms1_ptr_ = std::make_shared<MsAlignWriter>("msalign_1_" + std::to_string(i) + ".ms");
    MsAlignWriterPtr ms2_ptr_ = std::make_shared<MsAlignWriter>("msalign_2_" + std::to_string(i) + ".ms");
  
    writer_vec_1.push_back(ms1_ptr_);
    writer_vec_2.push_back(ms2_ptr_);
  }

  while (ms_group_ptr != nullptr) {
    while(pool_ptr->getQueueSize() >= thread_num * 2){
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    MatchEnvPtrVec prec_envs;

    RawMsGroupReader::obtainPrecEnvs(ms_group_ptr, prec_envs, 
                                     env_para_ptr_->prec_deconv_interval_, 
                                     env_para_ptr_->max_charge_); 

    EnvParaPtr env_ptr = deconv_ptr->getEnvPara();
    DpParaPtr dp_ptr = deconv_ptr->getDpPara();

    std::shared_ptr<EnvPara> env_ptr_new = std::make_shared<EnvPara>(*env_ptr);
    std::shared_ptr<DpPara> dp_ptr_new = std::make_shared<DpPara>(*dp_ptr);

    DeconvOneSpPtr deconv_ptr_new = std::make_shared<DeconvOneSp>(env_ptr_new, dp_ptr_new);
    
    MsAlignWriterPtr writer_ptr_1 = writer_vec_1[0];
    MsAlignWriterPtr writer_ptr_2 = writer_vec_2[0];

    writer_vec_1.erase(writer_vec_1.begin());
    writer_vec_2.erase(writer_vec_2.begin());

    pool_ptr->Enqueue(geneTask(ms_group_ptr, deconv_ptr_new, prec_envs, writer_ptr_1, writer_ptr_2, deconv_instance_ptr, &count1_ptr, &count2_ptr, total_scan_num, ms1_ptr_map_ptr, ms2_ptr_map_ptr));
    
    writer_vec_1.push_back(writer_ptr_1);
    writer_vec_2.push_back(writer_ptr_2);

    ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  }
  pool_ptr->ShutDown();

  writeMsalign(ms1_writer_ptr, ms1_ptr_map);
  writeMsalign(ms2_writer_ptr, ms2_ptr_map);

    //auto proc_end = std::chrono::high_resolution_clock::now();
    //auto proc_duration = std::chrono::duration_cast<std::chrono::microseconds>(proc_end-proc_start);
    //std::cout << std::endl << "Process " << proc_duration.count() << std::endl;

    //auto start = std::chrono::high_resolution_clock::now();
   
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    //std::cout << std::endl << "Read file " << duration.count() << std::endl;
  
  
  std::cout << std::endl;

  }
          
}; // namespace toppic
