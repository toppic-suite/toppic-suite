#include <boost/algorithm/string.hpp>

#include "log4cxx/logger.h"

#include "msalign_reader.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("MsAlignReader"));

MsAlignReader::MsAlignReader (const char *spectrum_file, 
                              ActivationPtrVec activation_list) {
  input_.open(spectrum_file, std::ios::in);
  activation_list_ = activation_list;
}

std::vector<std::string> MsAlignReader::readOneSpectrum() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    boost::algorithm::trim(line);
    if (line ==  "BEGIN IONS") {
      line_list.push_back(line);
    }
    else if (line == "END IONS") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    }
    else if (line == "") {
      continue;
    }
    else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

void MsAlignReader::readNext() {
  deconv_ms_ptr_ = DeconvMsPtr(nullptr);
  spectrum_str_ = readOneSpectrum();
  if (spectrum_str_.size() == 0) {
    input_.close();
    return;
  }
  std::vector<std::string> strs;
  int id = -1;
  std::string scans;
  std::string activation;
  std::string title;
  double prec_mass = -1;
  int prec_charge = -1;
  for (unsigned int i = 1; i < spectrum_str_.size() - 1; i++) {
    std::string letter = spectrum_str_[i].substr(0,1);
    if (letter >= "A" && letter <= "Z") {
      boost::split(strs, spectrum_str_[i], boost::is_any_of("="));
      if (strs[0] == "ID") {
        id = atoi(strs[1].c_str());
      }
      else if (strs[0] == "SCANS") {
        scans = strs[1];
      }
      else if (strs[0] == "ACTIVATION") {
        activation = strs[1];
      }
      else if (strs[0] == "TITLE") {
        title = strs[1];
      }
      else if (strs[0] == "PRECURSOR_MASS") {
        prec_mass = atof(strs[1].c_str());
      }
      else if (strs[0] == "PRECURSOR_CHARGE") {
        prec_charge = atoi(strs[1].c_str());
      }
    }
  }
  if (id < 0 || prec_charge < 0 || prec_mass < 0) {
    LOG4CXX_ERROR(logger, 
                  "Input file format error: sp id " << id << " prec_chrg "
                  << prec_charge << " prec mass " << prec_mass);
    std::exit(1);
  }

  MsHeaderPtr header_ptr(new MsHeader(prec_charge));
  header_ptr->setId(id);
  if (title != "") {
    std::stringstream ss;
    ss << "sp_" << id;
    header_ptr->setTitle(ss.str());
  } else {
    header_ptr->setTitle(title);
  }
  if (scans != "") {
    header_ptr->setScans(scans);
  }
  else {
    header_ptr->setScans("");
  }

  if (activation != "") {
    ActivationPtr activation_ptr = 
        getActivationPtrByName(activation_list_, activation);
    header_ptr->setActivationPtr(activation_ptr);
  }
  header_ptr->setPrecMonoMz(prec_mass /prec_charge
                       + MassConstant::getProtonMass());

  std::vector<DeconvPeakPtr> peak_ptr_list;
  int idx = 0;
  for (unsigned int i = 1; i < spectrum_str_.size() - 1; i++) {
    std::string letter = spectrum_str_[i].substr(0,1);
    if (letter >= "0" && letter <= "9") {
      boost::split(strs, spectrum_str_[i], boost::is_any_of("="));
      double mass = atof(strs[0].c_str());
      double inte = atof(strs[1].c_str());
      int charge = atoi(strs[2].c_str());
      DeconvPeakPtr peak_ptr(new DeconvPeak(idx, mass, inte, charge));
      peak_ptr_list.push_back(peak_ptr);
      idx++;
    }
  }
  deconv_ms_ptr_ = DeconvMsPtr(new Ms<DeconvPeakPtr>(header_ptr, peak_ptr_list));
  current_++;
}

DeconvMsPtr MsAlignReader::getNextMs() {
  readNext();
  return deconv_ms_ptr_;
}
    
void MsAlignReader::close() {
  input_.close();
}

int countSpNum(const char *spectrum_file, ActivationPtrVec activation_list) {
  MsAlignReader reader(spectrum_file, activation_list);
  int cnt = 0;
  DeconvMsPtr deconv_ms_ptr;
  while ((deconv_ms_ptr = reader.getNextMs()).get() != nullptr) {
    cnt++;
  }
  reader.close();
  return cnt;
}

}

