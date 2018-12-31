//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include "common/util/version.hpp"
#include "common/util/time_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"
#include "spec/msalign_writer.hpp"
#include "deconv/env/env_base.hpp"
#include "deconv/env/match_env_util.hpp"
#include "deconv/env/match_env_writer.hpp"
#include "deconv/deconv/deconv_process.hpp"

namespace toppic {

void DeconvProcess::copyParameters(EnvParaPtr env_para_ptr) {
  env_para_ptr->max_charge_ = para_ptr_->max_charge_;
  env_para_ptr->max_mass_ = para_ptr_->max_mass_;
  env_para_ptr->setTolerance(para_ptr_->tolerance_);
  env_para_ptr->ms_two_sn_ratio_ = para_ptr_->ms_two_sn_ratio_;
  env_para_ptr->ms_one_sn_ratio_ = para_ptr_->ms_one_sn_ratio_;
  env_para_ptr->keep_unused_peaks_ = para_ptr_->keep_unused_peaks_;
  env_para_ptr->output_multiple_mass_ = para_ptr_->output_multiple_mass_;
  env_para_ptr->prec_deconv_interval_ = para_ptr_->prec_window_;
  env_para_ptr->do_final_filtering_ = para_ptr_->do_final_filtering_;
}

std::string DeconvProcess::getParameterStr(DeconvParaPtr para_ptr, const std::string & prefix) {
  std::stringstream output;
  output << prefix << "TopFD " << version_number << std::endl;
  // TIME_STAMP_STR is replaced later
  output << prefix << "Timestamp: " << time_util::TIME_STAMP_STR << std::endl;
  output << prefix << "********************** Parameters **********************" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Input file: " << para_ptr->data_file_name_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Data type: " << "centroid" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Maximum charge: " << para_ptr->max_charge_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Maximum monoisotopic mass: " << para_ptr->max_mass_ << " Dalton" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Error tolerance: " << para_ptr->tolerance_ << " m/z" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "MS1 signal/noise ratio: " << para_ptr->ms_one_sn_ratio_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "MS/MS signal/noise ratio: " << para_ptr->ms_two_sn_ratio_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Precursor window size: " << para_ptr->prec_window_ << " m/z" << std::endl;
  //output << prefix << std::setw(40) << std::left 
  //    << "Do final filtering: " << para_ptr->do_final_filtering_ << std::endl;
  output << prefix << "********************** Parameters **********************" << std::endl;
  return output.str();
}

std::string DeconvProcess::updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num) {
  std::string percentage = str_util::toString(scan * 100 / total_scan_num);
  std::string msg = "Processing spectrum " + header_ptr->getTitle() + "...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

void DeconvProcess::process() {
  EnvBase::initBase(para_ptr_->resource_dir_);
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
  DpParaPtr dp_para_ptr = std::make_shared<DpPara>();
  copyParameters(env_para_ptr);
  std::string para_str = getParameterStr(para_ptr_);
  time_util::addTimeStamp(para_str);
  std::cout << para_str;

  std::string file_name = para_ptr_->getDataFileName();
  // writer
  std::string ms1_msalign_name, ms2_msalign_name;
  ms1_msalign_name = file_util::basename(file_name) + "_ms1.msalign";
  ms2_msalign_name = file_util::basename(file_name) + "_ms2.msalign";

  if (file_util::exists(file_util::basename(file_name) + "_ms1.env")) {
    file_util::delFile(file_util::basename(file_name) + "_ms1.env"); 
  }

  if (file_util::exists(file_util::basename(file_name) + "_ms2.env")) {
    file_util::delFile(file_util::basename(file_name) + "_ms2.env"); 
  }

  std::ofstream ms1_msalign_of(ms1_msalign_name, std::ofstream::out);
  std::ofstream ms2_msalign_of(ms2_msalign_name, std::ofstream::out);
  ms1_msalign_of.precision(16);
  ms2_msalign_of.precision(16);
  para_str = getParameterStr(para_ptr_, "#");
  time_util::addTimeStamp(para_str);
  ms1_msalign_of << para_str;
  ms2_msalign_of << para_str;

  DeconvOneSpPtr deconv_ptr = std::make_shared<DeconvOneSp>(env_para_ptr, dp_para_ptr);

  RawMsReaderPtr reader_ptr = std::make_shared<RawMsReader>(file_name);
  processSp(deconv_ptr, reader_ptr, ms1_msalign_of, ms2_msalign_of);

  ms1_msalign_of.close();
  ms2_msalign_of.close();
}

void DeconvProcess::processSp(DeconvOneSpPtr deconv_ptr, RawMsReaderPtr reader_ptr,
                              std::ofstream & ms1_msalign_of, std::ofstream & ms2_msalign_of) {
  // reader_ptr
  int total_scan_num = reader_ptr->getInputSpNum();
  RawMsPtr ms_ptr;
  int count1 = 0;
  int count2 = 0;
  ms_ptr = reader_ptr->getNextMs(para_ptr_->prec_window_, para_ptr_->max_charge_); 
  while (ms_ptr!= nullptr) {
    //auto proc_start = std::chrono::high_resolution_clock::now();
    PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
    LOG_DEBUG("peak list size " << peak_list.size());

    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    LOG_DEBUG("ms level " << header_ptr->getMsLevel() );

    // int scan_num_ = header_ptr->getFirstScanNum();
    std::string msg = updateMsg(header_ptr, count1 + count2 + 1, total_scan_num);
    std::cout << "\r" << msg << std::flush;
    LOG_DEBUG("set data....");
    deconv_ptr->setMsLevel(header_ptr->getMsLevel());
    if (header_ptr->getMsLevel() == 1) {
      deconv_ptr->setData(peak_list);
      deconv_ptr->run();
      MatchEnvPtrVec result_envs = deconv_ptr->getResult();
      LOG_DEBUG("result num " << result_envs.size());
      DeconvMsPtr ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);
      msalign_writer::write(ms1_msalign_of, ms_ptr);
      if (para_ptr_->output_match_env_) {
        match_env_writer::write(file_util::basename(para_ptr_->getDataFileName()) + "_ms1.env", header_ptr, result_envs);
      }
      count1++;
    } else {
      if (para_ptr_->missing_level_one_) {
        header_ptr->setPrecCharge(EnvPara::getDefaultMaxCharge());
        double prec_mz = EnvPara::getDefaultMaxMass()/EnvPara::getDefaultMaxCharge();
        header_ptr->setPrecMonoMz(prec_mz);
        header_ptr->setPrecSpMz(prec_mz);
        deconv_ptr->setData(peak_list, EnvPara::getDefaultMaxMass(),
                            EnvPara::getDefaultMaxCharge());
      } else {
        double max_frag_mass = header_ptr->getPrecMonoMass();
        if (max_frag_mass == 0.0) {
          max_frag_mass = header_ptr->getPrecSpMass();
        }
        deconv_ptr->setData(peak_list, max_frag_mass, header_ptr->getPrecCharge());
      }
      deconv_ptr->run();
      MatchEnvPtrVec result_envs = deconv_ptr->getResult();
      DeconvMsPtr ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);
      msalign_writer::write(ms2_msalign_of, ms_ptr);
      if (para_ptr_->output_match_env_) {
        match_env_writer::write(file_util::basename(para_ptr_->getDataFileName()) + "_ms2.env", header_ptr, result_envs);
      }
      count2++;
    }

    //auto proc_end = std::chrono::high_resolution_clock::now();
    //auto proc_duration = std::chrono::duration_cast<std::chrono::microseconds>(proc_end-proc_start);
    //std::cout << std::endl << "Process " << proc_duration.count() << std::endl;

    //auto start = std::chrono::high_resolution_clock::now();
    ms_ptr = reader_ptr->getNextMs(para_ptr_->prec_window_, para_ptr_->max_charge_);
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    //std::cout << std::endl << "Read file " << duration.count() << std::endl;
  }
  std::cout << std::endl << "Deconvolution finished." << std::endl;
}

}  // namespace toppic
