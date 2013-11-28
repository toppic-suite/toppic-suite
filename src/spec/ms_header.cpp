#include <sstream>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>

#include "log4cxx/logger.h"

#include "peak.hpp"
#include "ms_header.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("MsHeader"));

MsHeader::MsHeader(int chrg) {
  prec_chrg_ = chrg;
}


MsHeader::MsHeader(int scan_num, int level, int chrg) {
  scans_.push_back(scan_num);
  level_ = level;
  prec_chrg_ = chrg;
}

double MsHeader::getPrecMonoMass() {
  if (prec_mono_mz_ < 0) {
    LOG4CXX_WARN(logger, "monoisotopic mass is not initialized");
    return 0.0;
  } 
  else {
    return compPeakMass(prec_mono_mz_, prec_chrg_);
  }
}

double MsHeader::getPrecSpMass() {
  if (prec_sp_mz_ < 0) {
    LOG4CXX_WARN(logger, "precursor spectrum mass is not initialized");
    return 0.0;
  } else {
    return compPeakMass(prec_sp_mz_, prec_chrg_);
  }
}

double MsHeader::getPrecMonoMassMinusWater() {
  if (prec_mono_mz_ < 0) {
    LOG4CXX_WARN(logger, "monoisotopic mass is not initialized");
    return 0.0;
  } else {
    return compPeakMass(prec_mono_mz_, prec_chrg_)
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
		tmp << "Precursor Charge = " << prec_chrg_ << "\n";
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
  boost::split(strs,s,boost::is_any_of(" "));
  for (unsigned int i = 0; i < strs.size(); i++) {
    scans_.push_back(atoi(strs[i].c_str()));
  }
}

}
