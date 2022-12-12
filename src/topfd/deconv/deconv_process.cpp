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
#include <mutex>
#include <vector>
#include <boost/filesystem.hpp>

#include "common/util/logger.hpp"
#include "common/util/version.hpp"
#include "common/util/time_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/spec/msalign_thread_merge.hpp"

#include "ms/env/env_base.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_writer.hpp"

#include "topfd/common/topfd_para.hpp"
#include "topfd/msreader/raw_ms_writer.hpp"
#include "topfd/envcnn/env_cnn.hpp" 
#include "topfd/deconv/deconv_process.hpp"

namespace toppic {

std::mutex count_lock;
std::mutex count_lock_missing_ms1;
int DeconvProcess::ms1_spec_num_ = 0;
int DeconvProcess::ms2_spec_num_ = 0;

DeconvProcess::DeconvProcess(TopfdParaPtr topfd_para_ptr, 
                             const std::string &spec_file_name, 
                             int frac_id, int t) {
  topfd_para_ptr_ = topfd_para_ptr;
  env_para_ptr_ = std::make_shared<EnvPara>(topfd_para_ptr);
  dp_para_ptr_ = std::make_shared<DpPara>();

  spec_file_name_ = spec_file_name;
  frac_id_ = frac_id; 
  thread_num_ = t;
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

void DeconvProcess::prepareFileFolder(std::string file_num) {
  std::string base_name = file_util::basename(spec_file_name_);
  /*
  ms1_env_name_ = base_name_ + "_ms1.env"; 
  ms2_env_name_ = base_name_ + "_ms2.env"; 
  if (file_util::exists(ms1_env_name_)) {
    file_util::delFile(ms1_env_name_); 
  }
  if (file_util::exists(ms2_env_name_)) {
    file_util::delFile(ms2_env_name_); 
  }
  */

  if (topfd_para_ptr_->isGeneHtmlFolder()){
    //json file names
    html_dir_ =  base_name + "_" + file_num + "html";
    ms1_json_dir_ = html_dir_ 
        + file_util::getFileSeparator() + "topfd" 
        + file_util::getFileSeparator() + "ms1_json";
    ms2_json_dir_ = html_dir_ 
        + file_util::getFileSeparator() + "topfd" 
        + file_util::getFileSeparator() + "ms2_json";
    if (!file_util::exists(html_dir_)) {
      file_util::createFolder(html_dir_);
    }
    if (!file_util::exists(ms1_json_dir_)) {
      file_util::createFolder(ms1_json_dir_);
    }
    if (!file_util::exists(ms2_json_dir_)) {
      file_util::createFolder(ms2_json_dir_);
    }
  }
}

void DeconvProcess::process() {
  // reader
  RawMsGroupFaimeReaderPtr reader_ptr 
    = std::make_shared<RawMsGroupFaimeReader>(spec_file_name_, 
                                              topfd_para_ptr_->isMissingLevelOne(),
                                              topfd_para_ptr_->getActivation(),
                                              env_para_ptr_->prec_deconv_interval_, 
                                              frac_id_);
  //check if it is centroid data
  std::vector<pwiz::msdata::DataProcessingPtr> dt_proc_ptr = reader_ptr->getReaderPtr()->getMsdPtr()->dataProcessingPtrs;
  
  bool is_profile_data = true;

  if (dt_proc_ptr.size() > 0) {
    for (size_t i = 0; i < dt_proc_ptr.size(); i++) {
      for (size_t j = 0; j < dt_proc_ptr[i]->processingMethods.size(); j++) {
        pwiz::msdata::ProcessingMethod proc_method = dt_proc_ptr[i]->processingMethods[j];
        pwiz::cv::CVID cvid_peak_picking;
        cvid_peak_picking = pwiz::cv::MS_peak_picking;
        if (proc_method.hasCVParam(cvid_peak_picking)) {
          is_profile_data = false;
        }
      }
    }
  }
  if (is_profile_data) {
    LOG_ERROR("You are not using a centroid data.");
    LOG_ERROR("TopFD supports centroid data only. Please try again with a different dataset.");
    exit(EXIT_FAILURE);
  }

  if (topfd_para_ptr_->isUseEnvCnn()) {
    env_cnn::initModel(topfd_para_ptr_->getResourceDir(), topfd_para_ptr_->getThreadNum());
  }
  if (topfd_para_ptr_->isMissingLevelOne()) {
    processSpMissingLevelOne(reader_ptr);
  }
  else {
    processSp(reader_ptr);
  }
}

void deconvMissingMsOne(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, 
                        MsAlignWriterPtrVec ms_writer_ptr_vec, 
                        SimpleThreadPoolPtr pool_ptr, 
                        bool gene_html_folder, 
                        std::string ms2_json_dir){
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();

  header_ptr->setPrecCharge(EnvPara::getDefaultMaxCharge());
  double prec_mz = EnvPara::getDefaultMaxMass()/EnvPara::getDefaultMaxCharge();
  header_ptr->setPrecMonoMz(prec_mz);
  header_ptr->setPrecSpMz(prec_mz);
  MatchEnvPtrVec result_envs; 

  EnvParaPtr env_para_ptr = deconv_ptr->getEnvParaPtr();
  if (peak_list.size() > 0) {
    deconv_ptr->setData(peak_list, env_para_ptr->getMaxMass(), 
                        env_para_ptr->getMaxCharge());
    deconv_ptr->run();
    result_envs = deconv_ptr->getResult();
  }

  count_lock.lock();
  DeconvProcess::ms2_spec_num_++;
  count_lock.unlock();

  DeconvMsPtr deconv_ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);

  boost::thread::id thread_id = boost::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);

  ms_writer_ptr_vec[writer_id]->write(deconv_ms_ptr);

  /*
  if (topfd_para_ptr_->output_match_env_) {
    match_env_writer::write(ms2_env_name_, header_ptr, result_envs);
  }
  */

  if (gene_html_folder) {
    std::string json_file_name = ms2_json_dir 
        + file_util::getFileSeparator() 
        + "spectrum" + std::to_string(header_ptr->getId()) + ".js";
    raw_ms_writer::write(json_file_name, ms_ptr, result_envs);    
  }
}

std::function<void()> geneTaskMissingMsOne(RawMsGroupPtr ms_group_ptr, 
                                           DeconvOneSpPtr deconv_ptr,
                                           MsAlignWriterPtrVec ms_writer_ptr_vec, 
                                           SimpleThreadPoolPtr pool_ptr, 
                                           bool gene_html_dir,
                                           std::string ms2_json_dir){
  return [ms_group_ptr, deconv_ptr, ms_writer_ptr_vec, pool_ptr, 
  gene_html_dir, ms2_json_dir]() {

    RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();                           
    for (size_t i = 0; i < ms_two_ptr_vec.size(); i++) {
      RawMsPtr ms_two_ptr = ms_two_ptr_vec[i];
      deconvMissingMsOne(ms_two_ptr, deconv_ptr, ms_writer_ptr_vec, 
                         pool_ptr, gene_html_dir, ms2_json_dir);
    }
  };
}

void DeconvProcess::processSpMissingLevelOne(RawMsGroupFaimeReaderPtr reader_ptr) {
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num_);  

  MsAlignWriterPtrVec ms_writer_ptr_vec;
  std::string ms_msalign_name = file_util::basename(spec_file_name_) + "_ms2.msalign";
  for (int i = 0; i < thread_num_; i++) { 
    MsAlignWriterPtr ms_ptr = 
        std::make_shared<MsAlignWriter>(ms_msalign_name + "_" + str_util::toString(i));
    ms_writer_ptr_vec.push_back(ms_ptr);
  }

  // reader_ptr
  int total_scan_num = reader_ptr->getInputSpNum();
  RawMsGroupPtr ms_group_ptr;
  //ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  ms_group_ptr = reader_ptr->getNextMsGroupPtrWithFaime();

  if (topfd_para_ptr_->isGeneHtmlFolder()) {
    prepareFileFolder("");
  }

  int count = 0;
  while (ms_group_ptr != nullptr) {
    while(pool_ptr->getQueueSize() >= thread_num_ * 2){
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }

    //Here we create new instances of Env Para and Dp Para 
    //to make sure that multiple threads do not share these parameter instances. 
    EnvParaPtr env_ptr_new = std::make_shared<EnvPara>(*env_para_ptr_.get());
    DpParaPtr dp_ptr_new = std::make_shared<DpPara>(dp_para_ptr_);
    DeconvOneSpPtr deconv_ptr = std::make_shared<DeconvOneSp>(env_ptr_new, dp_ptr_new);

    pool_ptr->Enqueue(geneTaskMissingMsOne(ms_group_ptr, deconv_ptr, 
                                           ms_writer_ptr_vec, pool_ptr, 
                                           topfd_para_ptr_->isGeneHtmlFolder(), 
                                           ms2_json_dir_));

    RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();
    for (size_t t = 0; t < ms_two_ptr_vec.size(); t++){
      RawMsPtr ms_two_ptr = ms_two_ptr_vec[t];
      std::string msg = updateMsg(ms_two_ptr->getMsHeaderPtr(), count + 1, total_scan_num);
      std::cout << "\r" << msg << std::flush;
      count += 1;
    }
    //ms_group_ptr = reader_ptr->getNextMsGroupPtr();
    ms_group_ptr = reader_ptr->getNextMsGroupPtrWithFaime();

  }
  pool_ptr->ShutDown();

  /*MsalignThreadMerge msalign_ms2_merge = MsalignThreadMerge(
      file_util::basename(spec_file_name_), "ms2_msalign", 
      thread_num_, "ms2.msalign", topfd_para_ptr_-> getParaStr("#"));*/

  MsalignThreadMergePtr ms2_merge_ptr 
      = std::make_shared<MsalignThreadMerge>(file_util::basename(spec_file_name_), "ms2.msalign", 
                                             thread_num_, "ms2.msalign", 
                                             topfd_para_ptr_-> getParaStr("#", "\t"));

  ms2_merge_ptr->process();

  std::cout << std::endl;
}

void deconvMsOne(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, 
                 MatchEnvPtrVec prec_envs, 
                 MsAlignWriterPtrVec ms1_writer_ptr_vec, 
                 SimpleThreadPoolPtr pool_ptr, 
                 bool gene_html_dir,
                 std::string ms1_json_dir) { 
  
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  LOG_DEBUG("peak list size " << peak_list.size());
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  LOG_DEBUG("ms level " << header_ptr->getMsLevel() );
  // int scan_num_ = header_ptr->getFirstScanNum();

  count_lock.lock();
  DeconvProcess::ms1_spec_num_++;
  count_lock.unlock();
  
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

  boost::thread::id thread_id = boost::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);

  ms1_writer_ptr_vec[writer_id]->write(deconv_ms_ptr);

  /*
  if (topfd_para_ptr_->output_match_env_) {
    match_env_writer::write(ms1_env_name_, header_ptr, prec_envs);
  }
  */

  // add back precursor peaks
  for (size_t i = 0; i < peak_list.size(); i++) {
    peak_list[i]->setIntensity(intensities[i]);
  }
  //write only when html folder argument is true
  if (gene_html_dir) {
    std::string json_file_name = ms1_json_dir 
        + file_util::getFileSeparator() 
        + "spectrum" + std::to_string(header_ptr->getId()) + ".js";
    raw_ms_writer::write(json_file_name, ms_ptr, prec_envs);    
  }
}

void deconvMsTwo(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, 
                 MsAlignWriterPtrVec ms2_writer_ptr_vec, 
                 SimpleThreadPoolPtr pool_ptr, 
                 bool gene_html_dir, 
                 std::string ms2_json_dir) { 

  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  LOG_DEBUG("peak list size " << peak_list.size());
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  LOG_DEBUG("ms level " << header_ptr->getMsLevel() );
  // int scan_num_ = header_ptr->getFirstScanNum();
  double max_frag_mass = header_ptr->getPrecMonoMass();
  if (max_frag_mass == 0.0) {
    max_frag_mass = header_ptr->getPrecSpMass();
  }

  count_lock.lock();
  DeconvProcess::ms2_spec_num_++;
  count_lock.unlock();
  
  MatchEnvPtrVec result_envs;

  if (peak_list.size() > 0) {
    deconv_ptr->setMsLevel(header_ptr->getMsLevel());
    deconv_ptr->setData(peak_list, max_frag_mass, header_ptr->getPrecCharge());
    deconv_ptr->run();
    result_envs = deconv_ptr->getResult();
  }
  DeconvMsPtr deconv_ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);

  boost::thread::id thread_id = boost::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);

  ms2_writer_ptr_vec[writer_id]->write(deconv_ms_ptr);
  
  /*
  if (topfd_para_ptr_->output_match_env_) {
    match_env_writer::write(ms2_env_name_, header_ptr, result_envs);
  }
  */

  if (gene_html_dir) {
    std::string json_file_name = ms2_json_dir 
        + file_util::getFileSeparator() 
        + "spectrum" + std::to_string(ms_ptr->getMsHeaderPtr()->getId()) + ".js";
    raw_ms_writer::write(json_file_name, ms_ptr, result_envs);    
  }

}

//DecovOne & Two
std::function<void()> geneTask(RawMsGroupPtr ms_group_ptr, 
                               DeconvOneSpPtr deconv_ptr, 
                               MsAlignWriterPtrVec ms1_writer_ptr_vec, 
                               MsAlignWriterPtrVec ms2_writer_ptr_vec, 
                               SimpleThreadPoolPtr pool_ptr, 
                               bool gene_html_dir, 
                               std::string ms1_json_dir, 
                               std::string ms2_json_dir){
  return [ms_group_ptr, deconv_ptr, ms1_writer_ptr_vec, 
  ms2_writer_ptr_vec, pool_ptr, gene_html_dir, ms1_json_dir, ms2_json_dir]() {

    RawMsPtr ms_one_ptr = ms_group_ptr->getMsOnePtr();                            
    RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();

    MatchEnvPtrVec prec_envs;
    EnvParaPtr env_para_ptr = deconv_ptr->getEnvParaPtr();

    RawMsGroupFaimeReader::obtainPrecEnvs(ms_group_ptr, prec_envs, 
                                          env_para_ptr->max_mass_,
                                          env_para_ptr->max_charge_); 

    deconvMsOne(ms_one_ptr, deconv_ptr, prec_envs, ms1_writer_ptr_vec, 
                pool_ptr, gene_html_dir, ms1_json_dir);

    for (size_t i = 0; i < ms_two_ptr_vec.size(); i++) {
      RawMsPtr ms_two_ptr = ms_two_ptr_vec[i];
      deconvMsTwo(ms_two_ptr, deconv_ptr, ms2_writer_ptr_vec, 
                  pool_ptr, gene_html_dir, ms2_json_dir);
    }

  };
}

void DeconvProcess::processSp(RawMsGroupFaimeReaderPtr reader_ptr) {
  int total_scan_num = reader_ptr->getInputSpNum();//this is spectrum count, not same as scan ID count

  RawMsGroupPtr ms_group_ptr;

  ms_group_ptr = reader_ptr->getNextMsGroupPtrWithFaime();

  int count = 0;
  int vec_idx = 0; //index of the vector containing MsAlignWriterPtr for current msgroup
  
  std::string ms1_msalign_name, ms2_msalign_name;
  
  ms1_msalign_name = file_util::basename(spec_file_name_) + "_0_ms1.msalign";
  ms2_msalign_name = file_util::basename(spec_file_name_) + "_0_ms2.msalign";

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num_);  
  
  if (ms_group_ptr == nullptr) {
    LOG_ERROR("No spectrum to read in mzML file!");
    return;
  }
  //if it is not a FAIME data, don't add a number to a file name
  if (ms_group_ptr->getMsOnePtr()->getMsHeaderPtr()->getVoltage() == -1) {
    ms1_msalign_name = file_util::basename(spec_file_name_) + "_ms1.msalign";
    ms2_msalign_name = file_util::basename(spec_file_name_) + "_ms2.msalign";
    prepareFileFolder("");
  }
  else {
    ms1_msalign_name = file_util::basename(spec_file_name_) + "_0_ms1.msalign";
    ms2_msalign_name = file_util::basename(spec_file_name_) + "_0_ms2.msalign";
    prepareFileFolder("0_");
    isFaims_ = true;
  }

  //don't add voltage to vector if it is non-FAIME 
  voltage_vec_.push_back(std::make_pair(ms_group_ptr->getMsOnePtr()->getMsHeaderPtr()->getVoltage(), 0));

  std::vector<MsAlignWriterPtrVec> all_file_ms1_writer_ptr_vec;
  std::vector<MsAlignWriterPtrVec> all_file_ms2_writer_ptr_vec;

  MsAlignWriterPtrVec single_file_ms1_writer_ptr_vec;
  MsAlignWriterPtrVec single_file_ms2_writer_ptr_vec;
  
  //generate vector that contains msalign writing information
  MsAlignWriterPtrVec ms1_writer_ptr_vec;
  MsAlignWriterPtrVec ms2_writer_ptr_vec;

  for (int i = 0; i < thread_num_; i++) { 
    MsAlignWriterPtr ms1_ptr 
        = std::make_shared<MsAlignWriter>(ms1_msalign_name + "_" + str_util::toString(i));
    MsAlignWriterPtr ms2_ptr 
        = std::make_shared<MsAlignWriter>(ms2_msalign_name + "_" + str_util::toString(i));
    single_file_ms1_writer_ptr_vec.push_back(ms1_ptr);
    single_file_ms2_writer_ptr_vec.push_back(ms2_ptr);
  }
  all_file_ms1_writer_ptr_vec.push_back(single_file_ms1_writer_ptr_vec);
  all_file_ms2_writer_ptr_vec.push_back(single_file_ms2_writer_ptr_vec);
  while (ms_group_ptr != nullptr) {
    while(pool_ptr->getQueueSize() >= thread_num_ * 2){
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    bool is_new_voltage = true;
    //Here we create new Env Para and Dp Para instances to make sure that
    //multiple threads do not share the same parameter instances. 
    EnvParaPtr env_ptr_new = std::make_shared<EnvPara>(*env_para_ptr_.get());
    DpParaPtr dp_ptr_new = std::make_shared<DpPara>(dp_para_ptr_);

    //deconv_process_ptr_ (DecovProcess instance) is needed because 
    //it has the information on the folder names, envelope file names
    //pool_ptr needed for getting each thread id    
    DeconvOneSpPtr deconv_ptr = std::make_shared<DeconvOneSp>(env_ptr_new, dp_ptr_new);

    if (ms_group_ptr != nullptr) {
      //check if the voltage from this msgroup is new or not to determine whether to create a new set of writer vectors
      double cur_voltage = ms_group_ptr->getMsOnePtr()->getMsHeaderPtr()->getVoltage();
      for (size_t i = 0; i < voltage_vec_.size(); i++) {
        if (voltage_vec_[i].first == ms_group_ptr->getMsOnePtr()->getMsHeaderPtr()->getVoltage()) {
          int spec_id = voltage_vec_[i].second;
          RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();
          for (size_t j = 0; j < ms_two_ptr_vec.size(); j++) {
            ms_two_ptr_vec[j]->getMsHeaderPtr()->setMsOneId(spec_id);
          }
          voltage_vec_[i].second++;
          is_new_voltage = false;
          vec_idx = i;
          break;
        }
      }
      if (is_new_voltage) {
        //new file needs to be created
        std::string ms1_name, ms2_name;
        std::string file_num = str_util::toString(voltage_vec_.size());
    
        ms1_name = file_util::basename(spec_file_name_) + "_" + file_num + "_ms1.msalign";
        ms2_name = file_util::basename(spec_file_name_) + "_" + file_num + "_ms2.msalign";
        prepareFileFolder(file_num + "_");

        voltage_vec_.push_back(std::make_pair(cur_voltage, 0));

        int spec_id = 0;

        RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();

        for (size_t j = 0; j < ms_two_ptr_vec.size(); j++) {
          ms_two_ptr_vec[j]->getMsHeaderPtr()->setMsOneId(spec_id);
        }
        voltage_vec_.back().second++;

        single_file_ms1_writer_ptr_vec.clear();
        single_file_ms2_writer_ptr_vec.clear();

        for (int i = 0; i < thread_num_; i++) { 
          MsAlignWriterPtr ms1_ptr 
              = std::make_shared<MsAlignWriter>(ms1_name + "_" + str_util::toString(i));
          MsAlignWriterPtr ms2_ptr 
              = std::make_shared<MsAlignWriter>(ms2_name + "_" + str_util::toString(i));
          single_file_ms1_writer_ptr_vec.push_back(ms1_ptr);
          single_file_ms2_writer_ptr_vec.push_back(ms2_ptr);
        }
        all_file_ms1_writer_ptr_vec.push_back(single_file_ms1_writer_ptr_vec);
        all_file_ms2_writer_ptr_vec.push_back(single_file_ms2_writer_ptr_vec);

        vec_idx = all_file_ms1_writer_ptr_vec.size() -1;
      }
    }
    ms1_writer_ptr_vec = all_file_ms1_writer_ptr_vec[vec_idx];
    ms2_writer_ptr_vec = all_file_ms2_writer_ptr_vec[vec_idx];

    pool_ptr->Enqueue(geneTask(ms_group_ptr, deconv_ptr, ms1_writer_ptr_vec, 
                               ms2_writer_ptr_vec, pool_ptr, topfd_para_ptr_->isGeneHtmlFolder(), 
                               ms1_json_dir_, ms2_json_dir_));

    //count is 1 scan from msalign1 + n scan from msalign2 vector
    int parsed_scan = 1 + static_cast<int>((ms_group_ptr->getMsTwoPtrVec()).size());
    std::string msg = updateMsg(ms_group_ptr->getMsOnePtr()->getMsHeaderPtr(), count + parsed_scan, total_scan_num);
    std::cout << "\r" << msg << std::flush;
    count += parsed_scan;
    ms_group_ptr = reader_ptr->getNextMsGroupPtrWithFaime();    
  }
  pool_ptr->ShutDown();

  //auto proc_end = std::chrono::high_resolution_clock::now();
  ///auto proc_duration = std::chrono::duration_cast<std::chrono::microseconds>(proc_end-proc_start);
  //std::cout << std::endl << "Process " << proc_duration.count() << std::endl;

  //auto start = std::chrono::high_resolution_clock::now();

  //auto end = std::chrono::high_resolution_clock::now();
  //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
  //std::cout << std::endl << "Read file " << duration.count() << std::endl;
  std::string ms1_file_name = "ms1.msalign";
  std::string ms2_file_name = "ms2.msalign";

  if (voltage_vec_[0].first == -1) {
    //non-FAIME dataset
    MsalignThreadMergePtr ms1_merge_ptr
            = std::make_shared<MsalignThreadMerge>(spec_file_name_, ms1_file_name, 
                                                  thread_num_, ms1_file_name, 
                                                  topfd_para_ptr_-> getParaStr("#", "\t"));

    MsalignThreadMergePtr ms2_merge_ptr 
          = std::make_shared<MsalignThreadMerge>(spec_file_name_, ms2_file_name, 
                                                thread_num_, ms2_file_name, 
                                                topfd_para_ptr_-> getParaStr("#", "\t"));
    ms1_merge_ptr->process();
    ms2_merge_ptr->process();
  }
  else {
    for (size_t j = 0; j < voltage_vec_.size(); j++) {
      ms1_file_name = str_util::toString(j) + "_ms1.msalign";
      ms2_file_name = str_util::toString(j) + "_ms2.msalign";

      MsalignThreadMergePtr ms1_merge_ptr
            = std::make_shared<MsalignThreadMerge>(spec_file_name_, ms1_file_name, 
                                                  thread_num_, ms1_file_name, 
                                                  topfd_para_ptr_-> getParaStr("#", "\t"));

      MsalignThreadMergePtr ms2_merge_ptr 
            = std::make_shared<MsalignThreadMerge>(spec_file_name_, ms2_file_name, 
                                                  thread_num_, ms2_file_name, 
                                                  topfd_para_ptr_-> getParaStr("#", "\t"));
      ms1_merge_ptr->process();
      ms2_merge_ptr->process();

      msalign_num_ = voltage_vec_.size();
    }
  }
  
  std::cout << std::endl;
}

}; // namespace toppic
