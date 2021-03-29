//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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


#ifndef TOPPIC_MS_FEATURE_SAVITZKY_GOLAY_HPP_
#define TOPPIC_MS_FEATURE_SAVITZKY_GOLAY_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace toppic {

class SavitzkyGolay {
 public:
  SavitzkyGolay(int point_num, int poly_order);

  std::vector<double> smooth(std::vector<double> &values);

 private:
  int point_num_;
  boost::numeric::ublas::matrix<double> coef_mat_;
};

typedef std::shared_ptr<SavitzkyGolay> SavitzkyGolayPtr;

}

#endif
