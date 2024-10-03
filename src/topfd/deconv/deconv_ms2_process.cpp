// Copyright (c) 2014 - 2023, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "topfd/deconv/deconv_ms2_process.hpp"

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/feature/spec_feature_reader.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/mzml/mzml_ms_json_writer.hpp"
#include "ms/spec/msalign_thread_merge.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "topfd/deconv/deconv_single_sp.hpp"

namespace toppic {

// add a namespace to avoid duplicated method names
namespace deconv_ms2_process {

std::string updateMsTwoMsg(MsHeaderPtr header_ptr, int scan_cnt,
                           int total_scan_num) {
  std::string percentage = str_util::toString(scan_cnt * 100 / total_scan_num);
  std::string msg = "Processing MS/MS spectrum scan " +
                    std::to_string(header_ptr->getFirstScanNum()) + " ...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

void deconvMsTwo(MzmlMsPtr ms_ptr, SpecFeaturePtrVec sp_feat_ptr_vec,
                 TopfdParaPtr topfd_para_ptr,
                 MsAlignWriterPtrVec ms2_writer_ptr_vec,
                 SimpleThreadPoolPtr pool_ptr) {
  // 1. Find max_mass and max_charge
  double max_mass = 0;
  int max_charge = 1;
  for (std::size_t i = 0; i < sp_feat_ptr_vec.size(); i++) {
    double mass = sp_feat_ptr_vec[i]->getPrecMonoMass();
    if (mass > max_mass) {
      max_mass = mass;
    }
    int charge = sp_feat_ptr_vec[i]->getPrecCharge();
    if (charge > max_charge) {
      max_charge = charge;
    }
  }
  double arg_max_mass = topfd_para_ptr->getMaxMass();
  if (arg_max_mass < max_mass) {
    max_mass = arg_max_mass;
  }
  int arg_max_charge = topfd_para_ptr->getMaxCharge();
  if (arg_max_charge < max_charge) {
    max_charge = arg_max_charge;
  }
  // 2. Deconv the whole spectrum with filtering
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  MatchEnvPtrVec deconv_envs;
  if (peak_list.size() > 0) {
    int ms_level = 2;
    DeconvSingleSpPtr deconv_ptr = std::make_shared<DeconvSingleSp>(
        topfd_para_ptr, peak_list, ms_level, max_mass, max_charge);
    deconv_envs = deconv_ptr->deconv();
  }
  // 3. Add precuror information
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  std::sort(sp_feat_ptr_vec.begin(), sp_feat_ptr_vec.end(),
            SpecFeature::cmpPrecInteDec);
  PrecursorPtrVec prec_ptr_vec;
  for (std::size_t i = 0; i < sp_feat_ptr_vec.size(); i++) {
    int prec_id = i;
    int feat_id = sp_feat_ptr_vec[i]->getFracFeatureId();
    double mono_mz = sp_feat_ptr_vec[i]->getPrecMonoMz();
    int charge = sp_feat_ptr_vec[i]->getPrecCharge();
    double inte = sp_feat_ptr_vec[i]->getPrecInte();
    PrecursorPtr prec_ptr =
        std::make_shared<Precursor>(prec_id, feat_id, mono_mz, charge, inte);
    prec_ptr_vec.push_back(prec_ptr);
  }
  header_ptr->setPrecPtrVec(prec_ptr_vec);

  // 4. Write to msalign file
  DeconvMsPtr deconv_ms_ptr =
      match_env_util::getDeconvMsPtr(header_ptr, deconv_envs);

  boost::thread::id thread_id = boost::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);
  ms2_writer_ptr_vec[writer_id]->writeMs(deconv_ms_ptr);

  // 5. write json file
  if (topfd_para_ptr->isGeneHtmlFolder()) {
    std::string json_file_name =
        topfd_para_ptr->getMs2JsonDir() + file_util::getFileSeparator() +
        "spectrum" + std::to_string(header_ptr->getSpecId()) + ".js";
    mzml_ms_json_writer::write(json_file_name, ms_ptr, deconv_envs);
  }
}

std::function<void()> geneMsTwoTask(MzmlMsPtr ms_ptr,
                                    SpecFeaturePtrVec feat_ptr_vec,
                                    TopfdParaPtr topfd_para_ptr,
                                    MsAlignWriterPtrVec ms2_writer_ptr_vec,
                                    SimpleThreadPoolPtr pool_ptr) {
  return
      [ms_ptr, feat_ptr_vec, topfd_para_ptr, ms2_writer_ptr_vec, pool_ptr]() {
        deconvMsTwo(ms_ptr, feat_ptr_vec, topfd_para_ptr, ms2_writer_ptr_vec,
                    pool_ptr);
      };
}

}  // namespace deconv_ms2_process

DeconvMs2Process::DeconvMs2Process(TopfdParaPtr topfd_para_ptr, 
                                   const std::string & output_filename_ext) {
  topfd_para_ptr_ = topfd_para_ptr;
  output_filename_ext_ = output_filename_ext;
}

void DeconvMs2Process::prepareFileFolder() {
  if (topfd_para_ptr_->isGeneHtmlFolder()) {
    // json file names
    std::string html_dir = topfd_para_ptr_->getHtmlDir();
    if (!file_util::exists(html_dir)) {
      file_util::createFolder(html_dir);
    }
    std::string ms2_json_dir = topfd_para_ptr_->getMs2JsonDir();
    if (!file_util::exists(ms2_json_dir)) {
      file_util::createFolder(ms2_json_dir);
    }
  }
}

void DeconvMs2Process::readSpecFeature(
    std::string feat_file_name, std::map<int, SpecFeaturePtrVec> &feat_map) {
  SpecFeatureReaderPtr sp_feat_reader =
      std::make_shared<SpecFeatureReader>(feat_file_name);
  SpecFeaturePtrVec sp_feat_ptr_vec = sp_feat_reader->readAllFeatures();
  sp_feat_reader = nullptr;
  std::map<int, SpecFeaturePtrVec>::iterator feat_it;
  SpecFeaturePtrVec empty_vec;
  for (std::size_t i = 0; i < sp_feat_ptr_vec.size(); i++) {
    SpecFeaturePtr feat_ptr = sp_feat_ptr_vec[i];
    int sp_id = feat_ptr->getSpecId();
    feat_it = feat_map.find(sp_id);
    if (feat_it == feat_map.end()) {
      // if not found, insert
      feat_it = feat_map.insert(
          feat_map.end(), std::pair<int, SpecFeaturePtrVec>(sp_id, empty_vec));
    }
    feat_it->second.push_back(feat_ptr);
  }
}

void DeconvMs2Process::process() {
  MzmlMsGroupReaderPtr reader_ptr = std::make_shared<MzmlMsGroupReader>(
      topfd_para_ptr_->getMzmlFileName(), topfd_para_ptr_->getPrecWindowWidth(),
      topfd_para_ptr_->getActivation(), topfd_para_ptr_->getFracId(),
      topfd_para_ptr_->isFaims(), topfd_para_ptr_->getFaimsVoltage(),
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
  std::string output_base_name = topfd_para_ptr_->getOutputBaseName();
  std::string ms2_msalign_name = output_base_name + "_" + output_filename_ext_; 
  MsAlignWriterPtrVec ms2_writer_ptr_vec;
  for (int i = 0; i < thread_num; i++) {
    MsAlignWriterPtr ms2_ptr = std::make_shared<MsAlignWriter>(
        ms2_msalign_name + "_" + str_util::toString(i));
    ms2_writer_ptr_vec.push_back(ms2_ptr);
  }
  // reader spectrum features
  std::map<int, SpecFeaturePtrVec> feat_map;
  SpecFeaturePtrVec missing_one_feat_list;
  if (topfd_para_ptr_->isMissingLevelOne()) {
    double max_mass = topfd_para_ptr_->getMaxMass();
    int max_charge = topfd_para_ptr_->getMaxCharge();
    double mono_mz = peak_util::compMz(max_mass, max_charge);
    SpecFeaturePtr feat_ptr =
        std::make_shared<SpecFeature>(mono_mz, max_charge);
    missing_one_feat_list.push_back(feat_ptr);
  } else {
    std::string feat_file_name = output_base_name + "_ms2.feature";
    readSpecFeature(feat_file_name, feat_map);
  }
  // counter for processed spectra
  int spec_cnt = 0;
  // total spectrum number
  int total_spec_num = topfd_para_ptr_->getMs2ScanNum();
  std::map<int, SpecFeaturePtrVec>::iterator feat_it;
  while (ms_group_ptr != nullptr) {
    MzmlMsPtrVec ms_ptr_vec = ms_group_ptr->getMsTwoPtrVec();
    for (std::size_t i = 0; i < ms_ptr_vec.size(); i++) {
      while (pool_ptr->getQueueSize() >= thread_num * 2) {
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
      }
      spec_cnt++;
      MzmlMsPtr ms_ptr = ms_ptr_vec[i];
      SpecFeaturePtrVec sp_feat_ptr_vec;
      if (topfd_para_ptr_->isMissingLevelOne()) {
        sp_feat_ptr_vec = missing_one_feat_list;
      } else {
        feat_it = feat_map.find(ms_ptr->getMsHeaderPtr()->getSpecId());
        if (feat_it != feat_map.end()) {
          SpecFeaturePtrVec feat_list = feat_it->second;
          std::sort(feat_list.begin(), feat_list.end(), SpecFeature::cmpPrecInteDec);
          sp_feat_ptr_vec.push_back(feat_list[0]);
          double first_inte = feat_list[0]->getPrecInte();
          for (std::size_t i = 1; i < feat_list.size(); i++) {
            if (feat_list[i]->getPrecInte() >=
                first_inte * topfd_para_ptr_->getPrecInteCutoffRatio()) {
              sp_feat_ptr_vec.push_back(feat_list[i]);
              LOG_DEBUG("Inte " << feat_list[i]->getPrecInte() <<  " first inte " << first_inte); 
            }
          }
          LOG_DEBUG("Spectrum " << ms_ptr->getMsHeaderPtr()->getFirstScanNum()
                    << " feature " << feat_list.size() << " filtered " 
                    << (feat_list.size() -sp_feat_ptr_vec.size()) << " features."); 
        }
      }
      pool_ptr->Enqueue(deconv_ms2_process::geneMsTwoTask(
          ms_ptr, sp_feat_ptr_vec, topfd_para_ptr_, ms2_writer_ptr_vec,
          pool_ptr));
      std::string msg = deconv_ms2_process::updateMsTwoMsg(
          ms_ptr->getMsHeaderPtr(), spec_cnt, total_spec_num);
      std::cout << "\r" << msg << std::flush;
    }
    ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  }
  pool_ptr->ShutDown();
  for (int i = 0; i < thread_num; i++) {
    ms2_writer_ptr_vec[i] = nullptr;
  }
  // Merge files
  std::string para_str = topfd_para_ptr_->getParaStr("#", "\t");
  MsalignThreadMergePtr ms2_merge_ptr = std::make_shared<MsalignThreadMerge>(
      output_filename_ext_, topfd_para_ptr_->getThreadNum(), output_filename_ext_,
      output_base_name, para_str);
  ms2_merge_ptr->process();

  // remove tempory files
  std::string ms2_prefix =
      file_util::absoluteName(output_base_name) + "_" + output_filename_ext_ + "_";
  std::replace(output_base_name.begin(), output_base_name.end(), '\\', '/');
  file_util::cleanPrefix(output_base_name, ms2_prefix);
  std::cout << std::endl;
}

};  // namespace toppic
