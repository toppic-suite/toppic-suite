#ifndef PROT_COMP_SHIFT_LOW_MEM_HPP_
#define PROT_COMP_SHIFT_LOW_MEM_HPP_

#include <memory>
#include <vector>

namespace prot {

class CompShiftLowMem {
 public:

  std::vector<std::vector<int>> findBestShift(const std::vector<int> &a,
                                              const std::vector<int> &b);

  std::vector<double> findBestShift(const std::vector<int> &a,
                                    const std::vector<int> &errors,
                                    const std::vector<int> &b,
                                    int total,int minimum_gap,
                                    double scale);

  std::vector<std::vector<int>> findBestShift(const std::vector<int> &a,
                                              const std::vector<int> &errors,
                                              const std::vector<int> &b,
                                              int total,int minimum_gap) ;

 private:
  std::vector<short> num_;
  std::vector<std::vector<int>> findBestShift(const std::vector<int> &a,
                                              const std::vector<int> &b,
                                              int total,int minimum_gap) ;
  int checkD(std::vector<std::vector<int>> &ans,int d,int cur_min,
             int total,int min_gap);
};

typedef std::shared_ptr<CompShiftLowMem> CompShiftLowMemPtr;
} /* namespace prot */

#endif /* COMP_SHIFT_LOW_MEM_HPP_ */
