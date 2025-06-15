// Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane
// University.
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
//

#include "topfd/deconv/deconv_ms1_process.hpp"

#include <cstddef>
#include <memory>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/env_base.hpp"  // Ensure EnvBase is included
#include "ms/env/match_env_refine.hpp"  // Include EnvUtil for compEnvDist
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/mzml/mzml_ms_json_writer.hpp"
#include "ms/mzml/mzml_ms_sql_writer.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/spec/msalign_thread_merge.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "topfd/deconv/deconv_prec_win.hpp"
#include "topfd/deconv/deconv_single_sp.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

// add a namespace to avoid duplicated method names
namespace deconv_ms1_process {

class MassGroup {
 public:
  int cnt_ = 0;
  int offset_ = 0;
  double score_ = 0.0;
  MatchEnvPtrVec envs_;
};

typedef std::shared_ptr<MassGroup> MassGroupPtr;
typedef std::vector<MassGroupPtr> MassGroupPtrVec;

bool cmpMassGroupCntDecScoreDec(const MassGroupPtr &a, const MassGroupPtr &b) {
  if (a->cnt_ != b->cnt_) {
    return a->cnt_ > b->cnt_;
  }
  return a->score_ > b->score_;
}

std::string updateMsOneMsg(MsHeaderPtr header_ptr, int scan_cnt,
                           int total_scan_num) {
  std::string percentage = str_util::toString(scan_cnt * 100 / total_scan_num);
  std::string msg = "Processing MS1 spectrum scan " +
                    std::to_string(header_ptr->getFirstScanNum()) + " ...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

MatchEnvPtrVec refineMatchEnvVec(MatchEnvPtrVec envs,
                                 TopfdParaPtr topfd_para_ptr) {
  MatchEnvPtrVec result_envs;
  std::sort(envs.begin(), envs.end(), MatchEnv::cmpEnvcnnScoreDec);
  std::vector<double> offsets = topfd_para_ptr->getPrecOffsets();
  double tolerance = topfd_para_ptr->getPrecMergeTolerance();
  while (envs.size() > 0) {
    MassGroupPtrVec mass_groups;
    int base_mass_group = -1;
    for (size_t i = 0; i < offsets.size(); i++) {
      MassGroupPtr group_ptr = std::make_shared<MassGroup>();
      group_ptr->offset_ = offsets[i];
      mass_groups.push_back(group_ptr);
      if (offsets[i] == 0.0) {
        base_mass_group = i;
      }
    }
    // find an envelope group with similar masses
    mass_groups[base_mass_group]->cnt_++;
    mass_groups[base_mass_group]->score_ = envs[0]->getEnvcnnScore();
    mass_groups[base_mass_group]->envs_.push_back(envs[0]);
    double base_mass = envs[0]->getTheoEnvPtr()->getMonoNeutralMass();
    envs.erase(envs.begin());
    for (size_t i = envs.size() - 1; i >= 0; i--) {
      MatchEnvPtr env_ptr = envs[i];
      double env_mass = env_ptr->getTheoEnvPtr()->getMonoNeutralMass();
      double mass_diff = env_mass - base_mass;
      for (size_t j = 0; j < offsets.size(); j++) {
        double error = std::abs(mass_diff - offsets[j]);
        if (error < tolerance) {
          mass_groups[j]->envs_.push_back(env_ptr);
          envs.erase(envs.begin() + i);
          mass_groups[j]->cnt_++;
          if (mass_groups[j]->score_ < env_ptr->getEnvcnnScore()) {
            mass_groups[j]->score_ = env_ptr->getEnvcnnScore();
          }
          break;
        }
      }
    }
    // find the best mass group
    std::sort(mass_groups.begin(), mass_groups.end(),
              cmpMassGroupCntDecScoreDec);
    double best_offset = mass_groups[0]->offset_;
    result_envs.insert(result_envs.end(), mass_groups[0]->envs_.begin(),
                       mass_groups[0]->envs_.end());
    for (size_t i = 1; i < mass_groups.size(); i++) {
      double offset = mass_groups[i]->offset_;
      double mass_diff = best_offset - offset;
      int int_diff = static_cast<int>(mass_diff); 
      MatchEnvPtrVec &envs = mass_groups[i]->envs_;
      for (size_t j = 0; j < envs.size(); j++) {
        ExpEnvPtr real_env = envs[j]->getExpEnvPtr();
        double cur_mz = real_env->getReferMz();
        int charge = real_env->getCharge();
        double new_mz = cur_mz + mass_diff / charge;
        // check if the mass is greater than the precursor mass
        double ref_mass = peak_util::compPeakNeutralMass(new_mz, charge);
        // get a reference distribution based on the reference mass
        EnvPtr refer_env = EnvBase::getEnvByRefMass(ref_mass);
        /* add two zeros at both sides of the envelope */
        EnvPtr ext_refer_env = refer_env->addZero(2);
        // convert the reference distribution to a theoretical distribution
        // based on the base mz and charge state
        EnvPtr theo_env = ext_refer_env->distrToTheoRef(new_mz, charge);
        double max_inte = theo_env->getReferInte();
        theo_env->changeIntensity(1.0 / max_inte);
        int max_back_peak_num = real_env->getReferIdx()+int_diff;
        int max_forw_peak_num = 
            real_env->getPeakNum() - real_env->getReferIdx() - 1 - int_diff;
        MatchEnvPtr new_env_ptr;
        if (max_back_peak_num >= 0 && max_forw_peak_num >= 0) {
          EnvPtr new_theo_env = theo_env->getSubEnv(max_back_peak_num, max_forw_peak_num);
          double dist;
          double ratio;
          match_env_refine::compEnvDist(real_env, new_theo_env, dist, ratio);
          new_theo_env->changeIntensity(ratio);
          envs[j]->setTheoEnvPtr(new_theo_env);
          real_env->changeReferIdx(int_diff);
        }
        result_envs.push_back(envs[j]);
      }
    }
  }
  return result_envs;
}

void deconvMsOne(MzmlMsGroupPtr ms_group_ptr, TopfdParaPtr topfd_para_ptr,
                 MsAlignWriterPtrVec ms1_writer_ptr_vec,
                 SimpleThreadPoolPtr pool_ptr) {
  // 1. Store peak intensity
  MzmlMsPtr ms_ptr = ms_group_ptr->getMsOnePtr();
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  std::vector<double> intensities;
  for (std::size_t i = 0; i < peak_list.size(); i++) {
    intensities.push_back(peak_list[i]->getIntensity());
  }
  double base_inte = baseline_util::getBaseLine(intensities);
  double min_ref_inte = base_inte * topfd_para_ptr->getMsOneSnRatio();

  // 2. Deconv envelopes in precursor windows and remove them
  MatchEnvPtrVec prec_envs = deconv_prec_win::deconvPrecWinForMsGroup(
      ms_group_ptr, topfd_para_ptr->getMaxMass(),
      topfd_para_ptr->getMaxCharge(), base_inte, min_ref_inte);

  // Obtain EnvCNN Score for envelopes
  onnx_env_cnn::computeEnvScores(peak_list, prec_envs);

  // remove precursor peaks
  for (std::size_t i = 0; i < prec_envs.size(); i++) {
    ExpEnvPtr env_ptr = prec_envs[i]->getExpEnvPtr();
    for (int p = 0; p < env_ptr->getPeakNum(); p++) {
      if (env_ptr->isExist(p)) {
        int idx = env_ptr->getPeakIdx(p);
        peak_list[idx]->setIntensity(0);
      }
    }
  }
  // 3. Deconv the whole spectrum with filtering
  // get base intensity and min_ref_intensity for sql writing
  MatchEnvPtrVec deconv_envs;
  if (peak_list.size() > 0) {
    int ms_level = 1;
    double max_mass = topfd_para_ptr->getMaxMass();
    int max_charge = topfd_para_ptr->getMaxCharge();
    DeconvSingleSpPtr deconv_ptr = std::make_shared<DeconvSingleSp>(
        topfd_para_ptr, peak_list, ms_level, max_mass, max_charge);
    deconv_envs = deconv_ptr->deconv();
  }
  // 4. Merge precursor envelopes and deconvolution envelopes
  MatchEnvPtrVec result_envs;
  result_envs.insert(result_envs.end(), prec_envs.begin(), prec_envs.end());
  result_envs.insert(result_envs.end(), deconv_envs.begin(), deconv_envs.end());
  LOG_DEBUG("result num " << prec_envs.size());
  result_envs = refineMatchEnvVec(result_envs, topfd_para_ptr);

  // 5. Write to msalign file
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  DeconvMsPtr deconv_ms_ptr =
      match_env_util::getDeconvMsPtr(header_ptr, result_envs);

  boost::thread::id thread_id = boost::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);
  ms1_writer_ptr_vec[writer_id]->writeMs(deconv_ms_ptr);

  // 6. write json file
  if (topfd_para_ptr->isGeneHtmlFolder()) {
    // add back precursor peaks
    for (std::size_t i = 0; i < peak_list.size(); i++) {
      peak_list[i]->setIntensity(intensities[i]);
    }
    std::string json_file_name =
        topfd_para_ptr->getMs1JsonDir() + file_util::getFileSeparator() +
        "spectrum" + std::to_string(header_ptr->getSpecId()) + ".js";
    mzml_ms_json_writer::write(json_file_name, ms_ptr, prec_envs);
  }
  // 7. write sqlite file
  if (topfd_para_ptr->isGeneSql()) {
    LOG_DEBUG("Update SQl")
    mzml_ms_sql_writer::writeMs1(topfd_para_ptr->getSqlDb(), ms_ptr,
                                 result_envs, base_inte, min_ref_inte);
  }
}

std::function<void()> geneTask(MzmlMsGroupPtr ms_group_ptr,
                               TopfdParaPtr topfd_para_ptr,
                               MsAlignWriterPtrVec ms1_writer_ptr_vec,
                               SimpleThreadPoolPtr pool_ptr) {
  return [ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr]() {
    deconvMsOne(ms_group_ptr, topfd_para_ptr, ms1_writer_ptr_vec, pool_ptr);
  };
}

}  // namespace deconv_ms1_process

DeconvMs1Process::DeconvMs1Process(TopfdParaPtr topfd_para_ptr) {
  topfd_para_ptr_ = topfd_para_ptr;
}

void DeconvMs1Process::prepareFileFolder() {
  if (topfd_para_ptr_->isGeneHtmlFolder()) {
    // json file names
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

void DeconvMs1Process::process() {
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
  std::string ms1_msalign_name = output_base_name + "_ms1.msalign";
  MsAlignWriterPtrVec ms1_writer_ptr_vec;
  for (int i = 0; i < thread_num; i++) {
    MsAlignWriterPtr ms1_ptr = std::make_shared<MsAlignWriter>(
        ms1_msalign_name + "_" + str_util::toString(i));
    ms1_writer_ptr_vec.push_back(ms1_ptr);
  }
  // counter for processed spectra
  int spec_cnt = 0;
  // total spectrum number
  int total_spec_num = topfd_para_ptr_->getMs1ScanNum();
  while (ms_group_ptr != nullptr) {
    while (pool_ptr->getQueueSize() >= thread_num * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(deconv_ms1_process::geneTask(
        ms_group_ptr, topfd_para_ptr_, ms1_writer_ptr_vec, pool_ptr));
    spec_cnt++;
    std::string msg = deconv_ms1_process::updateMsOneMsg(
        ms_group_ptr->getMsOnePtr()->getMsHeaderPtr(), spec_cnt,
        total_spec_num);
    std::cout << "\r" << msg << std::flush;
    ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  }
  pool_ptr->ShutDown();
  for (int i = 0; i < thread_num; i++) {
    ms1_writer_ptr_vec[i] = nullptr;
  }
  // Merge files
  std::string file_name_ext = "ms1.msalign";
  std::string para_str = topfd_para_ptr_->getParaStr("#", "\t");
  MsalignThreadMergePtr ms1_merge_ptr = std::make_shared<MsalignThreadMerge>(
      file_name_ext, topfd_para_ptr_->getThreadNum(), file_name_ext,
      output_base_name, para_str);
  ms1_merge_ptr->process();

  // remove tempory files
  std::string ms1_prefix =
      file_util::absoluteName(output_base_name) + "_ms1.msalign_";
  std::replace(output_base_name.begin(), output_base_name.end(), '\\', '/');
  file_util::cleanPrefix(output_base_name, ms1_prefix);
  std::cout << std::endl;
}

};  // namespace toppic
