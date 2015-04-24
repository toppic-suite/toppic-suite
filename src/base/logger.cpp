#include "base/logger.hpp"

namespace prot {

int log_level = 5;

std::ofstream WebLog::log_;
std::string WebLog::log_file_name_;
double WebLog::ratio_;

void WebLog::init(std::string log_file_name, bool use_gf) {
  log_file_name_ = log_file_name;
  if (log_file_name_.length() > 0) {
    log_.open(log_file_name, std::ios::out | std::ios::app);
  }
  if (use_gf) {
	  ratio_ = 1.0;  
  } else {
    ratio_ = 3.5;
  }
}

void WebLog::percentLog(double p) {
  if (log_.is_open()) {
    log_ << p * ratio_ << std::endl;
  }
}

void WebLog::close() {
  if (log_.is_open()) {
    log_ << 1 << std::endl;
    log_.close();
  }
}

}
