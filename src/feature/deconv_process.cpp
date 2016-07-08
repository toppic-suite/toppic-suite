#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "feature/deconv_process.hpp"
#include "feature/msalign_writer.hpp"

namespace prot {

void DeconvProcess::copyParameters(FeatureMngPtr mng_ptr) {
  mng_ptr->max_charge_ = para_ptr_->max_charge_;
  mng_ptr->max_mass_ = para_ptr_->max_mass_;
  mng_ptr->setTolerance(para_ptr_->tolerance);
  mng_ptr->sn_ratio_ = para_ptr_->sn_ratio_;
  mng_ptr->keep_unused_peaks_ = para_ptr_->keep_unused_peaks_;
  mng_ptr->output_multiple_mass_ = para_ptr_->output_multiple_mass_;
  mng_ptr->prec_deconv_interval_ = para_ptr_->prec_window_;
}

void DeconvProcess::printParameter(FeatureMngPtr mng_ptr) {
  std::cout << "TopFD " << std::endl;
  std::cout << "********* parameters begin **********" << std::endl;
  std::cout << "output file format:    " <<  para_ptr_->output_type_ << std::endl;
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
  FeatureMngPtr mng_ptr(new FeatureMng());
  copyParameters(mng_ptr);
  printParameter(mng_ptr);

  std::string file_name = para_ptr_->getDataFileName();
  // writer
  std::string base_name = FileUtil::basename(file_name) + ".msalign";
  std::ofstream of(base_name, std::ofstream::out);

  DeconvOneSpPtr deconv_ptr(new DeconvOneSp(mng_ptr));

  FeatureMsReaderPtr reader_ptr(new FeatureMsReader(file_name));
  processSp(deconv_ptr, reader_ptr, of);
}

void DeconvProcess::processSp(DeconvOneSpPtr deconv_ptr, FeatureMsReaderPtr reader_ptr, 
                              std::ofstream &of) {
  // reader_ptr
  int total_scan_num = reader_ptr->getInputSpNum();
  RawMsPtr ms_ptr;
  int count = 0;
  while ((ms_ptr = reader_ptr->getNextMs(para_ptr_->prec_window_)) != nullptr) {
    PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
    if (peak_list.size() == 0) {
      continue;
    }
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    if (header_ptr->getMsLevel() == 1 &&  para_ptr_->ms_level_ != 1) {
      continue;
    }
    int scan_num_ = header_ptr->getFirstScanNum();
    updateMsg(header_ptr, scan_num_, total_scan_num);
    LOG_DEBUG("set data....");
    if (header_ptr->getMsLevel() == 1) {
      deconv_ptr->setData(peak_list);
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
        double max_frag_mass = 0;
        for (size_t i = 0; i < peak_list.size(); i++) {
          if (peak_list[i]->getPosition() > max_frag_mass) {
            max_frag_mass = peak_list[i]->getPosition();
          }
        }
        
        deconv_ptr->setData(peak_list, max_frag_mass, header_ptr->getPrecCharge());
      }
    }
    deconv_ptr->run();
    MatchEnvPtrVec result_envs = deconv_ptr->getResult();
    MsalignWriter::writeText(of, result_envs, header_ptr, count);
    count++;
  }
  std::cout <<"Deconvolution finished." << std::endl;
}

}
