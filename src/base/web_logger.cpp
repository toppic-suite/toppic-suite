#include "base/web_logger.hpp"

namespace prot {

std::ofstream WebLog::log_;
std::string WebLog::log_file_name_;

double WebLog::total_time_;
double WebLog::function_start_time_;

void WebLog::init(std::string log_file_name, bool use_gf, bool localization, int ptm_num) {
  log_file_name_ = log_file_name;
  if (log_file_name_.length() > 0) {
    log_.open(log_file_name, std::ios::out | std::ios::app);
  }
  function_start_time_ = 0;
  total_time_ = ZeroPtmFilterTime() + ZeroPtmSearchTime();
  if (ptm_num >= 1) {
    total_time_ += OnePtmFilterTime();
    total_time_ += OnePtmSearchTime();
  }
  if (ptm_num >= 2) {
    total_time_ += DiagFilterTime();
    total_time_ += TwoPtmSearchTime();
  }
  if (use_gf) {
    total_time_ += GfEvalueTime();  
  } else {
    total_time_ += TableEvalueTime();
  }
  if (localization) {
    total_time_ += LocalizationTime();
  }
  total_time_ += OutPutTime();
}

void WebLog::close() {
  if (log_.is_open()) {
    log_ << 1 << std::endl;
    log_.close();
  }
}

void WebLog::percentLog(int spectrum_index, int spectrum_num, int block_index, 
                        int block_num, double func_time) {
  double func_percentage = (double)block_index/block_num 
      + (double)spectrum_index/spectrum_num/block_num;
  double time = function_start_time_ + func_time * func_percentage;
  double total_percentage = time / total_time_;
  if (log_.is_open()) {
    log_ << total_percentage << std::endl;
  }
}

void WebLog::percentLog(int spectrum_index, int spectrum_num, double func_time) {
  return percentLog(spectrum_index, spectrum_num, 0, 1, func_time);
}

void WebLog::completeFunction(double func_time) {
  function_start_time_ += func_time;
}

}
