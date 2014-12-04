#include "base/logger.hpp"

namespace prot {

int log_level = 5;

std::ofstream WebLog::log_;
std::string WebLog::log_file_name_;

void WebLog::init(std::string log_file_name) {
  log_file_name_ = log_file_name;
  if (log_file_name_.length() > 0) {
    log_.open(log_file_name, std::ios::out | std::ios::app);
  }
}

void WebLog::percent_log(double p) {
  if (log_.is_open()) {
    log_ << p << std::endl;
  }
}

void WebLog::close() {
  if (log_.is_open()) {
    log_ << 1 << std::endl;
    log_.close();
  }
}

}
