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


#ifndef TOPPIC_TOPLIB_LIB_NUMERIC_HPP_
#define TOPPIC_TOPLIB_LIB_NUMERIC_HPP_

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace toppic {

namespace toplib {

// The arithmetic of the Python libraries TopLib's ms_library_building.py
// runs on, reproduced so that the C++ port computes the same numbers.
namespace lib_numeric {

// numpy's sum of a float64 array (np.sum, ndarray.sum): the pairwise sum
// with eight accumulators on blocks of up to 128 values.
double numpySum(const double* a, size_t n);

inline double numpySum(const std::vector<double>& a) {
  return numpySum(a.data(), a.size());
}

// Python's built-in sum() of a list of floats (Python >= 3.12): Neumaier's
// compensated summation.
double pythonSum(const std::vector<double>& a);

// np.linalg.norm of a float64 vector: the square root of OpenBLAS's ddot,
// whose x86-64 AVX-512 kernel (SkylakeX, numpy's bundled OpenBLAS) sums the
// products in 32 lanes of fused multiply-adds, then 16, then one by one.
double numpyNorm(const std::vector<double>& a);

// A number of a table cell as pandas' read_csv parses it (the C parser's
// precise_xstrtod: up to 17 significant digits accumulated in a double, then
// scaled by an exact power of ten, which is not always correctly rounded;
// "inf" and "infinity" with an optional sign, in any case). False if the
// text is not a number.
bool pandasToDouble(const std::string& text, double& value);

// Python's random module (the Mersenne Twister) for the calls TopLib makes:
// random.seed(n) with a non-negative integer and random.sample(range(n), k)
// with a population of at most 21 items.
class PyRandom {
 public:
  explicit PyRandom(uint64_t seed);

  // random.sample(range(n), k), n <= 21
  std::vector<int> sample(int n, int k);

 private:
  uint32_t genrandUint32();

  uint32_t getrandbits(int k) { return genrandUint32() >> (32 - k); }

  // random._randbelow_with_getrandbits
  int randbelow(int n);

  void initGenrand(uint32_t s);

  void initByArray(const std::vector<uint32_t>& key);

  static constexpr int kN = 624;
  static constexpr int kM = 397;
  std::array<uint32_t, kN> mt_;
  int index_ = kN;
};

}  // namespace lib_numeric

}  // namespace toplib

}  // namespace toppic

#endif
