#include <cmath>

#include "base/logger.hpp"
#include "base/activation_base.hpp"
#include "base/string_util.hpp"
#include "spec/msalign_reader.hpp"

namespace prot {

MsAlignReader::MsAlignReader(const std::string &file_name, 
                             int group_spec_num, ActivationPtr act_ptr): 
    file_name_(file_name),
    group_spec_num_(group_spec_num),
    activation_ptr_(act_ptr) {
      input_.open(file_name.c_str(), std::ios::in);
    }

std::vector<std::string> MsAlignReader::readOneSpectrum() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    line = StringUtil::trim(line);
    if (line ==  "BEGIN IONS") {
      line_list.push_back(line);
    } else if (line == "END IONS") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } else if (line == "") {
      continue;
    } else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

void MsAlignReader::readNext() {
  deconv_ms_ptr_ = DeconvMsPtr(nullptr);
  spectrum_str_vec_ = readOneSpectrum();
  if (spectrum_str_vec_.size() == 0) {
    input_.close();
    return;
  }
  int id = -1;
  int prec_id = 0;
  std::string scans;
  double retention_time = -1;
  std::string activation;
  std::string title;
  int level = 2;
  double prec_mass = -1;
  int prec_charge = -1;
  double prec_inte = -1;
  std::vector<std::string> strs;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0,1);
    if (letter >= "A" && letter <= "Z") {
      strs = StringUtil::split(spectrum_str_vec_[i], '=');
      if (strs[0] == "ID") {
        id = std::stoi(strs[1]);
      }

      if (strs[0] == "PRECURSOR_ID") {
        prec_id = std::stoi(strs[1]);
      } else if (strs[0] == "SCANS") {
        scans = strs[1];
      } else if (strs[0] == "RETENTION_TIME") {
        retention_time = std::stod(strs[1]);        
      } else if (strs[0] == "ACTIVATION") {
        activation = strs[1];
      } else if (strs[0] == "TITLE") {
        title = strs[1];
      } else if (strs[0] == "LEVEL") {
        level = std::stoi(strs[1]);
      } else if (strs[0] == "PRECURSOR_MASS") {
        prec_mass = std::stod(strs[1]);
      } else if (strs[0] == "PRECURSOR_CHARGE") {
        prec_charge = std::stoi(strs[1]);
      } else if (strs[0] == "PRECURSOR_INTENSITY") {
        prec_inte = std::stod(strs[1]);
      }
    }
  }
  if (id < 0 || prec_charge < 0 || prec_mass < 0) {
    LOG_ERROR("Input file format error: sp id " << id << " prec_chrg "
              << prec_charge << " prec mass " << prec_mass);
    std::exit(1);
  }

  MsHeaderPtr header_ptr(new MsHeader());
  header_ptr->setFileName(file_name_);
  header_ptr->setId(id);
  header_ptr->setPrecId(prec_id);
  if (scans != "") {
    header_ptr->setScans(scans);
  } else {
    header_ptr->setScans("");
  }
  header_ptr->setRetentionTime(retention_time);
  //LOG_DEBUG("retention time " << retention_time);

  if (title != "") {
    std::stringstream ss;
    ss << "sp_" << id;
    header_ptr->setTitle(ss.str());
  } else {
    header_ptr->setTitle(title);
  }

  if (activation_ptr_ != nullptr) {
    header_ptr->setActivationPtr(activation_ptr_);
  } else if (activation != "") {
    ActivationPtr activation_ptr = 
        ActivationBase::getActivationPtrByName(activation);
    header_ptr->setActivationPtr(activation_ptr);
  }
  header_ptr->setMsLevel(level);

  header_ptr->setPrecMonoMz(prec_mass /prec_charge
                            + MassConstant::getProtonMass());
  header_ptr->setPrecCharge(prec_charge);
  
  header_ptr->setPrecInte(prec_inte);

  std::vector<DeconvPeakPtr> peak_ptr_list;
  int idx = 0;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0,1);
    if (letter >= "0" && letter <= "9") {
      strs = StringUtil::split(spectrum_str_vec_[i], '\t');
      double mass = std::stod(strs[0]);
      double inte = std::stod(strs[1]);
      int charge = std::stoi(strs[2]);
      DeconvPeakPtr peak_ptr(new DeconvPeak(idx, mass, inte, charge));
      peak_ptr_list.push_back(peak_ptr);
      idx++;
    }
  }
  deconv_ms_ptr_ 
      = DeconvMsPtr(new Ms<DeconvPeakPtr>(header_ptr, peak_ptr_list));

  current_++;
}

DeconvMsPtr MsAlignReader::getNextMs() {
  readNext();
  return deconv_ms_ptr_;
}

SpectrumSetPtr MsAlignReader::getNextSpectrumSet(SpParaPtr sp_para_ptr) {
  DeconvMsPtrVec deconv_ms_ptr_vec; 
  for (int i = 0; i < group_spec_num_; i++) {
    readNext();
    if (deconv_ms_ptr_ == nullptr) {
      return SpectrumSetPtr(nullptr);
    }
    deconv_ms_ptr_vec.push_back(deconv_ms_ptr_);
  }
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  //LOG_DEBUG("prec mass " << prec_mono_mass);
  int count = 1;
  for (int i = 1; i < group_spec_num_; i++) {
    double new_mass = deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getPrecMonoMass();
    if (std::abs(prec_mono_mass - new_mass) < 0.5) {
      prec_mono_mass = (prec_mono_mass * count + new_mass)/ (count+1);
      count++;
    }
  }
  //LOG_DEBUG("prec mass result " << prec_mono_mass);
  return SpectrumSetPtr(new SpectrumSet(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass));
}

void MsAlignReader::close() {
  input_.close();
}

}
