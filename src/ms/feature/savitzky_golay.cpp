//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "ms/feature/matrix_inverse.hpp"
#include "ms/feature/savitzky_golay.hpp"

namespace toppic {

namespace ublas = boost::numeric::ublas;

void getFilterMatrix(int poly_order, int filter_len, ublas::matrix<double> &mat) {
  // filter_len is always an odd number
  int m = (filter_len - 1) /2;
  ublas::matrix<double> a(filter_len, poly_order + 1);
  for (int i = -m; i <= m; i++) {
    for (int j = 0; j <= poly_order; j++) {
      double di = i;
      double dj = j;
      a(i + m,j) = std::pow(di,dj);
    }
  }
  ublas::matrix<double> a_t = ublas::trans(a); 
  ublas::matrix<double> f = ublas::prod(a_t, a);
  ublas::matrix<double> f_inv;

  matrix_inverse::InvertMatrix(f, f_inv);

  ublas::matrix<double> a_f_inv = ublas::prod(a, f_inv);
  mat = ublas::prod(a_f_inv, a_t);
}

bool isEmpty(std::vector<double> &values) {
  for (size_t i = 0; i < values.size(); i++) {
    if (values[i] > 0) {
      return false;
    }
  }
  return true;
}

SavitzkyGolay::SavitzkyGolay(int point_num, int poly_order) {
  if (point_num < 3 || point_num %2 == 0) {
    LOG_ERROR("Point number is too small or is not an odd number!");
    exit(EXIT_FAILURE);
  }
  point_num_ = point_num;
  getFilterMatrix(poly_order, point_num, coef_mat_); 
} 


std::vector<double> SavitzkyGolay::smooth(std::vector<double> &values) {
  // No need to smooth if all values are 0
  if (isEmpty(values)) {
    return values;
  }
  int m = (point_num_ - 1)/2;
  int cnt = values.size();
  std::vector<double> result(cnt, 0.0);

  for (int i = 0; i <= m; i++) {
    double dot_prod = 0;
    for (int j = 0; j < point_num_; j++) {
      dot_prod += coef_mat_(i, j) * values[j];
    }
    result[i] = dot_prod;
  }

  for (int i = m + 1; i < cnt - m - 1; i++) {
    double dot_prod = 0.0;
    for (int j = 0; j < point_num_; j++) {
      dot_prod += coef_mat_(m, j) * values[i-m + j];
    }
    result[i] = dot_prod;
  }

  for (int i = 0; i <= m; i++) {
    double dot_prod = 0;
    for (int j = 0; j < point_num_; j++) {
      dot_prod += coef_mat_(m+ i, j) * values[cnt - point_num_ +j];
    }
    result[cnt - m - 1 + i] = dot_prod;
  }

  return result;
}

} /* namespace */

