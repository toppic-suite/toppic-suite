//
// Created by abbash on 8/26/22.
//

#ifndef TOPPIC_UTILITY_FUNCTIONS_HPP
#define TOPPIC_UTILITY_FUNCTIONS_HPP

#include <vector>
namespace toppic {
namespace utility_functions {
  double pearsonr(std::vector<double> X, std::vector<double> Y);
  std::vector<double> findLocalMinima(std::vector<double> arr);
  std::vector<double> findLocalMaxima(std::vector<double> arr);
}
}


#endif //TOPPIC_UTILITY_FUNCTIONS_HPP
