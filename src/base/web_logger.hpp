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
  static double ZeroPtmFilterTime() {return 1.0;}
  static double ZeroPtmSearchTime() {return 1.0;}
  static double OnePtmFilterTime() {return 1.0;}
  static double OnePtmSearchTime() {return 2.0;}
  static double DiagFilterTime() {return 6.0;}
  static double TwoPtmSearchTime() {return 6.0;}
  static double TableEvalueTime() {return 2.0;}
  static double GfEvalueTime() {return 76.0;}
  static double SelectingTime() {return 1.0;}
  static double LocalizationTime() {return 3.0;}
  static double OutPutTime() {return 1.0;}

  static void init(std::string log_file_name, bool use_gf, bool localization, int ptm_num);
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
