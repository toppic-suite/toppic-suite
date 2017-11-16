//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include <string>
#include <iomanip>

#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "base/version.hpp"
#include "feature/deconv_process.hpp"
#include "feature/match_env.hpp"

namespace prot {

void DeconvProcess::copyParameters(FeatureMngPtr mng_ptr) {
  mng_ptr->max_charge_ = para_ptr_->max_charge_;
  mng_ptr->max_mass_ = para_ptr_->max_mass_;
  mng_ptr->setTolerance(para_ptr_->tolerance_);
  mng_ptr->sn_ratio_ = para_ptr_->sn_ratio_;
  mng_ptr->keep_unused_peaks_ = para_ptr_->keep_unused_peaks_;
  mng_ptr->output_multiple_mass_ = para_ptr_->output_multiple_mass_;
  mng_ptr->prec_deconv_interval_ = para_ptr_->prec_window_;
}

void DeconvProcess::outputParameter(std::ostream &output, DeconvParaPtr para_ptr, const std::string & prefix) {
  time_t cur_time = std::time(0);
  output << prefix << "TopFD " << version_number << std::endl;
  output << prefix << "Timestamp: " << asctime(localtime(&cur_time)); 
  output << prefix << "********************** Parameters **********************" << std::endl;
  output << prefix << std::setw(40) << std::left << "Input file: " << para_ptr->data_file_name_ << std::endl;
  output << prefix << std::setw(40) << std::left << "Data type: " << "centroided" << std::endl;
  output << prefix << std::setw(40) << std::left << "Maximum charge: " << para_ptr->max_charge_ << std::endl;
  output << prefix << std::setw(40) << std::left << "Maximum monoisotopic mass: " << para_ptr->max_mass_ << " Dalton" << std::endl;
  output << prefix << std::setw(40) << std::left << "Error tolerance: " << para_ptr->tolerance_ << " m/z" << std::endl;
  output << prefix << std::setw(40) << std::left << "Signal/noise ratio: " << para_ptr->sn_ratio_ << std::endl;
  output << prefix << std::setw(40) << std::left << "Precursor window size: " << para_ptr->prec_window_ << " m/z" << std::endl;
  output << prefix << "********************** Parameters **********************" << std::endl;
}

std::string DeconvProcess::updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num) {
  std::string percentage = std::to_string(scan * 100 / total_scan_num);
  std::string msg = "Processing spectrum " + header_ptr->getTitle() + "...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

void DeconvProcess::process() {
  FeatureMngPtr mng_ptr = std::make_shared<FeatureMng>(para_ptr_->exec_dir_);
  copyParameters(mng_ptr);
  outputParameter(std::cout, para_ptr_);

  std::string file_name = para_ptr_->getDataFileName();
  // writer
  std::string ms1_name, ms2_name;
  ms1_name = FileUtil::basename(file_name) + "_ms1.msalign";
  ms2_name = FileUtil::basename(file_name) + "_ms2.msalign";

  std::ofstream of1(ms1_name, std::ofstream::out);
  std::ofstream of2(ms2_name, std::ofstream::out);
  of1.precision(16);
  of2.precision(16);
  outputParameter(of1, para_ptr_, "#");
  outputParameter(of2, para_ptr_, "#");

  DeconvOneSpPtr deconv_ptr = std::make_shared<DeconvOneSp>(mng_ptr);

  FeatureMsReaderPtr reader_ptr = std::make_shared<FeatureMsReader>(file_name);
  processSp(deconv_ptr, reader_ptr, of1, of2);
  of1.close();
  of2.close();
}

void DeconvProcess::processSp(DeconvOneSpPtr deconv_ptr, FeatureMsReaderPtr reader_ptr,
                              std::ofstream &of1, std::ofstream &of2) {
  // reader_ptr
  int total_scan_num = reader_ptr->getInputSpNum();
  RawMsPtr ms_ptr;
  int count1 = 0;
  int count2 = 0;
  while ((ms_ptr = reader_ptr->getNextMs(para_ptr_->prec_window_)) != nullptr) {
    PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
    LOG_DEBUG("peak list size " << peak_list.size());

    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    LOG_DEBUG("ms level " << header_ptr->getMsLevel() );

    // int scan_num_ = header_ptr->getFirstScanNum();
    std::string msg = updateMsg(header_ptr, count1 + count2 + 1, total_scan_num);
    std::cout << "\r" << msg;
    LOG_DEBUG("set data....");
    deconv_ptr->setMsLevel(header_ptr->getMsLevel());
    if (header_ptr->getMsLevel() == 1) {
      deconv_ptr->setData(peak_list);
      deconv_ptr->run();
      MatchEnvPtrVec result_envs = deconv_ptr->getResult();
      MsalignWriter::write(of1, result_envs, header_ptr);
      count1++;
    } else {
      if (para_ptr_->missing_level_one_) {
        header_ptr->setPrecCharge(FeatureMng::getDefaultMaxCharge());
        double prec_mz = FeatureMng::getDefaultMaxMass()/FeatureMng::getDefaultMaxCharge();
        header_ptr->setPrecMonoMz(prec_mz);
        header_ptr->setPrecSpMz(prec_mz);
        deconv_ptr->setData(peak_list, FeatureMng::getDefaultMaxMass(),
                            FeatureMng::getDefaultMaxCharge());
      } else {
        double max_frag_mass = header_ptr->getPrecMonoMass();
        if (max_frag_mass == 0.0) {
          max_frag_mass = header_ptr->getPrecSpMass();
        }
        deconv_ptr->setData(peak_list, max_frag_mass, header_ptr->getPrecCharge());
      }
      deconv_ptr->run();
      MatchEnvPtrVec result_envs = deconv_ptr->getResult();
      MsalignWriter::write(of2, result_envs, header_ptr);
      count2++;
    }
  }
  std::cout << std::endl << "Deconvolution finished." << std::endl;
}

}  // namespace prot
