//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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

#include "topfd/ecscore/env_set/xic_util.hpp"

namespace toppic {

namespace xic_util {

std::vector<int> findLocalMinima(std::vector<double> &arr) {
  int n = arr.size();
  std::vector<int> minima;
  for (int i = 1; i < n - 1; i++) {
    if ((arr[i - 1] > arr[i]) and (arr[i] < arr[i + 1])) {
      if (i - 2 > 0)
        if (arr[i - 2] <= arr[i])
          continue;
      if (i + 2 < n)
        if (arr[i + 2] <= arr[i])
          continue;
      minima.push_back(i);
    }
  }
  return minima;
}

std::vector<int> findLocalMaxima(std::vector<double> &arr) {
  int n = arr.size();
  std::vector<int> maxima;
  for (int i = 1; i < n - 1; i++)
    if ((arr[i - 1] < arr[i]) and (arr[i] > arr[i + 1]))
      maxima.push_back(i);
  if (arr[n - 1] > arr[n - 2]) maxima.push_back(n - 1);
  return maxima;
}

}

}



