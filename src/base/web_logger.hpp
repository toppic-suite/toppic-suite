/*
 * author  Xiaowen Liu
 * date    2015-7-11
 */
#ifndef PROT_BASE_WEB_LOGGER_HPP_
#define PROT_BASE_WEB_LOGGER_HPP_

#include <string>
#include <iostream>
#include <fstream>

namespace prot {

class WebLog {
 public:
  static double ZeroPtmTime() {return 3.0;}
  static double DiagFilterTime() {return 5.0;}
  static double OnePtmFilterTime() {return 1.0;}
  static double PtmTime() {return 17.0;}
  static double TableEvalueTime() {return 2.0;}
  static double GfEvalueTime() {return 73.0;}

  static void init(std::string log_file_name, bool use_gf, int ptm_num);
  static void close();

  static void percentLog(int spectrum_index, int spectrum_num,
                         int block_index, int block_num, double func_time);

  static void percentLog(int spectrum_index, int spectrum_num, double func_time);

  static void completeFunction(double func_time);

 private:
  static std::string log_file_name_;
  static std::ofstream log_;

  static double total_time_;
  static double function_start_time_;
};

}
#endif
