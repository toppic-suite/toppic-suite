//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef TOPPIC_MS_UTIL_SAVITZKY_GOLAY_HPP_
#define TOPPIC_MS_UTIL_SAVITZKY_GOLAY_HPP_

#include <memory>
#include <vector>

// Boost uBLAS still derives its iterators from the C++17-deprecated
// std::iterator; silence that third-party warning so our build stays clean.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/numeric/ublas/matrix.hpp>
#pragma GCC diagnostic pop

namespace toppic {

class SavitzkyGolay {
 public:
  SavitzkyGolay(int point_num, int poly_order);

  std::vector<double> smooth(const std::vector<double> &values);

 private:
  int point_num_;
  boost::numeric::ublas::matrix<double> coef_mat_;
};

using SavitzkyGolayPtr = std::shared_ptr<SavitzkyGolay>;

}  // namespace toppic

#endif
