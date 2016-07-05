#include "base/logger.hpp"
#include "base/activation_base.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/common/cv.hpp"
#include "spec/raw_ms_reader.hpp" 

namespace prot {

RawMsReader::RawMsReader(std::string &file_name):
    file_name_(file_name),
    input_sp_id_(0),
    output_sp_id_(0) { 
      msd_ptr_ = MSDataFilePtr(new pwiz::msdata::MSDataFile(file_name, &readers_));
      spec_list_ptr_ =  msd_ptr_->run.spectrumListPtr;
      input_sp_num_ = spec_list_ptr_->size();
    }

void RawMsReader::readNext() {
  peak_list_.clear();
  header_ptr_ = nullptr;

  if (input_sp_id_ >= input_sp_num_) {
    return;
  }

  pwiz::msdata::SpectrumPtr cur_spec_ptr = nullptr;
  //Read m/z and intensity values from the spectra
  bool get_binary_data = true;	
  while (cur_spec_ptr == nullptr) {
    cur_spec_ptr = spec_list_ptr_->spectrum(input_sp_id_, get_binary_data);
    input_sp_id_++;
    if (input_sp_id_ >= input_sp_num_) {
      LOG_ERROR("Only " << input_sp_num_  << " spectra in the input data!");
      return;
    }
  }
  pwiz::msdata::SpectrumInfo spec_info(*cur_spec_ptr);
  for (size_t i = 0; i < spec_info.data.size(); i++) {
    PeakPtr peak_ptr (new Peak(spec_info.data[i].mz, spec_info.data[i].intensity));
    peak_list_.push_back(peak_ptr);
  }
  int ms_level = spec_info.msLevel; 
  if (ms_level == 2) {
    double prec_mz = spec_info.precursors[0].mz;
    if (prec_mz < 0) {
      prec_mz = 0;
    }
    int prec_charge = (int)spec_info.precursors[0].charge;
    if (prec_charge  < 0) {
      prec_charge = 1;
    }
    header_ptr_ = MsHeaderPtr(new MsHeader());
    header_ptr_->setScan(spec_info.scanNumber);
    header_ptr_->setMsLevel(ms_level);
    header_ptr_->setPrecCharge(prec_charge);
    header_ptr_->setFileName(file_name_);
    header_ptr_->setTitle("Scan_" + std::to_string(spec_info.scanNumber));
    // here is average mz 
    header_ptr_->setPrecSpMz(prec_mz);
    header_ptr_->setRetentionTime(spec_info.retentionTime);
    pwiz::cv::CVID ac_id = spec_info.massAnalyzerType;
    std::string ac_name;
    if (ac_id == pwiz::cv::MS_CID) {
      ac_name = "CID"; 
    }
    else if (ac_id == pwiz::cv::MS_HCD) {
      ac_name = "HCD";
    }
    else if (ac_id == pwiz::cv::MS_ETD) {
      ac_name = "ETD";
    }
    ActivationPtr activation_ptr = ActivationBase::getActivationPtrByName(ac_name); 
    header_ptr_->setActivationPtr(activation_ptr);
  } 
  else {
    header_ptr_ = MsHeaderPtr(new MsHeader());
    header_ptr_->setScan(spec_info.scanNumber);
    header_ptr_->setMsLevel(ms_level);
    header_ptr_->setPrecCharge(0);
    header_ptr_->setFileName(file_name_);
    header_ptr_->setTitle("Scan_" + std::to_string(spec_info.scanNumber));
    header_ptr_->setRetentionTime(spec_info.retentionTime);
  }
}

}
