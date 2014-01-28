/*
 * comp_shift_low_mem.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#ifndef PROT_COMP_SHIFT_LOW_MEM_HPP_
#define PROT_COMP_SHIFT_LOW_MEM_HPP_

#include <memory>
#include <vector>

namespace prot {

class CompShiftLowMem {
 public:
  CompShiftLowMem();
  const static int max_len_ = 1000;
  std::vector<std::vector<int>> findBestShift(std::vector<int> a,
                                              std::vector<int> b);
  std::vector<double> findBestShift(std::vector<int> a,
                                    std::vector<int> errors,
                                    std::vector<int> b,
                                    int total,int minimum_gap,
                                    double scale);
  std::vector<std::vector<int>> findBestShift(std::vector<int> a,
                                              std::vector<int> errors,
                                              std::vector<int> b,
                                              int total,int minimum_gap) ;
 protected:
  std::vector<short> getNum(){return num_;}

 private:
  std::vector<short> num_;
  std::vector<std::vector<int>> findBestShift(std::vector<int> a,
                                              std::vector<int> b,
                                              int total,int minimum_gap) ;
  int checkD(std::vector<std::vector<int>> &ans,int d,int cur_min,
             int total,int min_gap);
};

typedef std::shared_ptr<CompShiftLowMem> CompShiftLowMemPtr;
} /* namespace prot */

#endif /* COMP_SHIFT_LOW_MEM_HPP_ */
