// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "common/emass/isotope_data.hpp"

namespace toppic {

namespace emass {

const std::vector<Isotope>& getIsotopes(Element element) {
  static const std::array<std::vector<Isotope>, kElementNum> isotopes = {{
      // 12C, 13C
      {{12.0, 0.9893}, {13.00335483534, 0.0107}},
      // 1H, 2H
      {{1.007825031898, 0.999885}, {2.014101777844, 0.000115}},
      // 14N, 15N
      {{14.00307400425, 0.99636}, {15.00010889827, 0.00364}},
      // 16O, 17O, 18O
      {{15.99491461926, 0.99757}, {16.99913175595, 0.00038},
       {17.99915961214, 0.00205}},
      // 32S, 33S, 34S, no 35S, 36S
      {{31.97207117354, 0.9499}, {32.97145890862, 0.0075},
       {33.967867011, 0.0425}, {0.0, 0.0}, {35.967080692, 0.0001}},
  }};
  return isotopes[element];
}

const std::string& getSymbol(Element element) {
  static const std::array<std::string, kElementNum> symbols = {
      "C", "H", "N", "O", "S"};
  return symbols[element];
}

}  // namespace emass

}  // namespace toppic
