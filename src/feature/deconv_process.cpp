#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "feature/deconv_process.hpp"
#include "feature/msalign_writer.hpp"

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

void DeconvProcess::printParameter(FeatureMngPtr mng_ptr) {
  std::cout << "TopFD 1.0.0 (" << __DATE__ << ")" << std::endl;
  std::cout << "********* parameters begin **********" << std::endl;
  std::cout << "output file format:    " << para_ptr_->output_type_ << std::endl;
  std::cout << "data type:             " << "centroided" << std::endl;
  std::cout << "maximum charge:        " << mng_ptr->max_charge_ << std::endl;
  std::cout << "maximum mass:          " << mng_ptr->max_mass_ << std::endl;
  std::cout << "m/z error tolerance:   " << mng_ptr->mz_tolerance_ << std::endl;
  std::cout << "sn ratio:              " << mng_ptr->sn_ratio_ << std::endl;
  std::cout << "keep unused peak:      " << mng_ptr->keep_unused_peaks_ << std::endl;
  std::cout << "output multiple mass:  " << mng_ptr->output_multiple_mass_ << std::endl;
  std::cout << "********* parameters end   **********" << std::endl;
}

void DeconvProcess::updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num) {
  std::string percentage = std::to_string(scan * 100 / total_scan_num);
  msg_ = "Processing spectrum " + header_ptr->getTitle() + "...";
  while (msg_.length() < 40) {
    msg_ += " ";
  }
  msg_ = msg_ + percentage + "% finished.";
}

void DeconvProcess::process() {
  FeatureMngPtr mng_ptr(new FeatureMng(para_ptr_->exec_dir_));
  copyParameters(mng_ptr);
  printParameter(mng_ptr);

  std::string file_name = para_ptr_->getDataFileName();
  // writer
  std::string ms1_name = FileUtil::basename(file_name) + ".ms1";
  std::string ms2_name = FileUtil::basename(file_name) + ".ms2";
  std::ofstream of1(ms1_name, std::ofstream::out);
  std::ofstream of2(ms2_name, std::ofstream::out);
  of1.precision(16);
  of2.precision(16);

  DeconvOneSpPtr deconv_ptr(new DeconvOneSp(mng_ptr));

  FeatureMsReaderPtr reader_ptr(new FeatureMsReader(file_name));
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
    LOG_DEBUG(" peak list size " << peak_list.size());
    if (peak_list.size() == 0) {
      continue;
    }
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    LOG_DEBUG("ms level " << header_ptr->getMsLevel() );
    if (header_ptr->getMsLevel() == 1 &&  para_ptr_->ms_level_ != 1) {
      continue;
    }
    //int scan_num_ = header_ptr->getFirstScanNum();
    updateMsg(header_ptr, count1 + count2 + 1, total_scan_num);
    std::cout << "\r" << msg_;
    LOG_DEBUG("set data....");
    if (header_ptr->getMsLevel() == 1) {
      deconv_ptr->setData(peak_list);
      deconv_ptr->run();
      MatchEnvPtrVec result_envs = deconv_ptr->getResult();
      MsalignWriter::writeText(of1, result_envs, header_ptr);
      count1++;
    }
    else {
      if (para_ptr_->missing_level_one_) {
        header_ptr->setPrecCharge(FeatureMng::getDefaultMaxCharge());
        double prec_mz = FeatureMng::getDefaultMaxMass()/FeatureMng::getDefaultMaxCharge();
        header_ptr->setPrecMonoMz(prec_mz);
        header_ptr->setPrecSpMz(prec_mz);
        deconv_ptr->setData(peak_list, FeatureMng::getDefaultMaxMass(), 
                            FeatureMng::getDefaultMaxCharge());
      }
      else {
        double max_frag_mass = header_ptr->getPrecMonoMass();
        if (max_frag_mass == 0.0) {
          max_frag_mass = header_ptr->getPrecSpMass();
        }
        deconv_ptr->setData(peak_list, max_frag_mass, header_ptr->getPrecCharge());
      }
      deconv_ptr->run();
      MatchEnvPtrVec result_envs = deconv_ptr->getResult();
      MsalignWriter::writeText(of2, result_envs, header_ptr);
      count2++;
    }
  }
  //std::cout <<"Deconvolution finished." << std::endl;
}

}
