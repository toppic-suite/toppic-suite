//
// Created by abbash on 8/26/22.
//

#include <cmath>
#include "utility_functions.hpp"

namespace toppic{
namespace utility_functions{
  double pearsonr(std::vector<double> X, std::vector<double> Y){
    int n = X.size();
    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double squareSum_X = 0, squareSum_Y = 0;
    for (int i = 0; i < n; i++){
      sum_X = sum_X + X[i];
      sum_Y = sum_Y + Y[i];
      sum_XY = sum_XY + X[i] * Y[i];
      squareSum_X = squareSum_X + X[i] * X[i];
      squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }
    double corr = (double)(n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
    return corr;
  }
}
}