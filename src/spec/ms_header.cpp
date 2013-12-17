#include <sstream>
#include <stdlib.h>

#include <base/logger.hpp>
#include <base/string_util.hpp>
#include "spec/peak.hpp"
#include "spec/ms_header.hpp"

namespace prot {


double MsHeader::getPrecMonoMass() {
  if (prec_mono_mz_ < 0) {
    LOG_WARN("monoisotopic mass is not initialized")
    return 0.0;
  } 
  else {
    return compPeakMass(prec_mono_mz_, prec_charge_);
  }
}

double MsHeader::getPrecSpMass() {
  if (prec_sp_mz_ < 0) {
    LOG_WARN("precursor spectrum mass is not initialized");
    return 0.0;
  } else {
    return compPeakMass(prec_sp_mz_, prec_charge_);
  }
}

double MsHeader::getPrecMonoMassMinusWater() {
  if (prec_mono_mz_ < 0) {
    LOG_WARN("monoisotopic mass is not initialized");
    return 0.0;
  } else {
    return compPeakMass(prec_mono_mz_, prec_charge_)
        - MassConstant::getWaterMass();
  }
}

std::string MsHeader::toString() {
  std::stringstream tmp;
		tmp << "MS HEADER\n";
		tmp << "==========\n";
		tmp << "Title = " << title_ << "\n";
		tmp << "Scan Number = " << scans_[0] << "\n";
		tmp << "MS Level = " << level_ << "\n";
		tmp << "Activation type = " << activation_ptr_ << "\n";
		tmp << "Precursor Sp Mz = " << prec_sp_mz_ << "\n";
		tmp << "Precursor Charge = " << prec_charge_ << "\n";
		tmp << "Precursro Mono Mz = " << prec_mono_mz_ << "\n";
		return tmp.str();
	}

std::string MsHeader::getScansString() {
  std::stringstream scan_list;
  scan_list << scans_[0];
  for (unsigned int i = 1; i < scans_.size(); i++) {
    scan_list <<  " " << scans_[i];
  }
  return scan_list.str();
}



void MsHeader::setScans(std::string s) {
  if (s == "") {
    scans_.clear();
    scans_.push_back(-1);
    return;
  }
  std::vector<std::string> strs;
  split(s,' ', strs); 
  for (unsigned int i = 0; i < strs.size(); i++) {
    scans_.push_back(atoi(strs[i].c_str()));
  }
}

}
