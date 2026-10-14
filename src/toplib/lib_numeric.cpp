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


#include "toplib/lib_numeric.hpp"

#include <charconv>
#include <cmath>
#include <cstdio>
#include <limits>

namespace toppic {

namespace toplib {

namespace lib_numeric {

// numpy's DOUBLE_pairwise_sum (umath/loops_utils.h)
double numpySum(const double* a, size_t n) {
  if (n < 8) {
    double res = 0.;
    for (size_t i = 0; i < n; i++) {
      res += a[i];
    }
    return res;
  }
  if (n <= 128) {
    double r[8];
    for (int j = 0; j < 8; j++) {
      r[j] = a[j];
    }
    size_t i = 8;
    for (; i < n - (n % 8); i += 8) {
      for (int j = 0; j < 8; j++) {
        r[j] += a[i + j];
      }
    }
    double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; i++) {
      res += a[i];
    }
    return res;
  }
  size_t n2 = n / 2;
  n2 -= n2 % 8;
  return numpySum(a, n2) + numpySum(a + n2, n - n2);
}

// CPython's builtin_sum_impl: the int start 0 plus the first float gives
// the first float, then Neumaier's summation over the rest
double pythonSum(const std::vector<double>& a) {
  if (a.empty()) {
    return 0;
  }
  double f = 0 + a[0];
  double c = 0.0;
  for (size_t i = 1; i < a.size(); i++) {
    double x = a[i];
    double t = f + x;
    if (std::fabs(f) >= std::fabs(x)) {
      c += (f - t) + x;
    } else {
      c += (x - t) + f;
    }
    f = t;
  }
  if (c != 0 && std::isfinite(c)) {
    f += c;
  }
  return f;
}

// OpenBLAS kernel/x86_64/ddot.c with ddot_microk_skylakex-2.c: blocks of 16
// values (32 per step in four 8-lane accumulators, the rest 16 per step in
// four 4-lane ones), the lanes added pairwise, then the remaining values
double numpyNorm(const std::vector<double>& x) {
  size_t n = x.size();
  double dot = 0.0;
  size_t n1 = n & ~static_cast<size_t>(15);
  if (n1 > 0) {
    double acc[4][8] = {};
    size_t i = 0;
    size_t n32 = n1 & ~static_cast<size_t>(31);
    for (; i < n32; i += 32) {
      for (int a = 0; a < 4; a++) {
        for (int l = 0; l < 8; l++) {
          double v = x[i + 8 * a + l];
          acc[a][l] = std::fma(v, v, acc[a][l]);
        }
      }
    }
    double y[4][4];
    for (int a = 0; a < 4; a++) {
      for (int l = 0; l < 4; l++) {
        y[a][l] = acc[a][l] + acc[a][l + 4];
      }
    }
    for (; i < n1; i += 16) {
      for (int a = 0; a < 4; a++) {
        for (int l = 0; l < 4; l++) {
          double v = x[i + 4 * a + l];
          y[a][l] = std::fma(v, v, y[a][l]);
        }
      }
    }
    double s[4];
    for (int l = 0; l < 4; l++) {
      s[l] = ((y[0][l] + y[1][l]) + y[2][l]) + y[3][l];
    }
    dot = (s[0] + s[2]) + (s[1] + s[3]);
  }
  for (size_t i = n1; i < n; i++) {
    dot = std::fma(x[i], x[i], dot);
  }
  return std::sqrt(dot);
}

std::vector<double> log2Normalize(const std::vector<double>& v) {
  std::vector<double> r(v.size());
  for (size_t i = 0; i < v.size(); i++) {
    r[i] = std::log2(v[i]);
  }
  double norm = numpyNorm(r);
  for (double& x : r) {
    x = x / norm;
  }
  return r;
}

// Python's float_repr_style 'short' (format_float_short, mode 'r'):
// std::to_chars gives the same shortest, nearest digits as dtoa's mode 0
std::string pythonRepr(double v) {
  if (std::isnan(v)) {
    return "nan";
  }
  if (std::isinf(v)) {
    return v > 0 ? "inf" : "-inf";
  }
  char buf[64];
  auto res = std::to_chars(buf, buf + sizeof(buf), v,
                           std::chars_format::scientific);
  std::string sci(buf, res.ptr);
  std::string sign;
  if (!sci.empty() && sci[0] == '-') {
    sign = "-";
    sci.erase(0, 1);
  }
  size_t e = sci.find('e');
  std::string digits;
  for (size_t i = 0; i < e; i++) {
    if (sci[i] != '.') {
      digits += sci[i];
    }
  }
  // the decimal point after decpt digits, as dtoa reports it
  int decpt = std::stoi(sci.substr(e + 1)) + 1;
  int n = static_cast<int>(digits.size());
  std::string out;
  if (decpt <= -4 || decpt > 16) {
    out = digits.substr(0, 1);
    if (n > 1) {
      out += "." + digits.substr(1);
    }
    int exp = decpt - 1;
    char eb[16];
    std::snprintf(eb, sizeof(eb), "e%+.02d", exp);
    out += eb;
  } else if (decpt <= 0) {
    out = "0." + std::string(-decpt, '0') + digits;
  } else if (decpt < n) {
    out = digits.substr(0, decpt) + "." + digits.substr(decpt);
  } else {
    out = digits + std::string(decpt - n, '0') + ".0";
  }
  return sign + out;
}

namespace {

bool isSpace(char c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' ||
         c == '\r';
}

bool isDigit(char c) { return c >= '0' && c <= '9'; }

// the powers of ten precise_xstrtod scales by
const std::array<double, 309>& powersOfTen() {
  static const std::array<double, 309> e = [] {
    std::array<double, 309> t{};
    for (int i = 0; i < 309; i++) {
      t[i] = std::stod("1e" + std::to_string(i));
    }
    return t;
  }();
  return e;
}

}  // namespace

// pandas/_libs/src/parser/tokenizer.c: precise_xstrtod (decimal '.', sci
// 'E', no thousands separator, trailing spaces skipped) as to_double calls
// it, and _try_double's spellings of infinity
bool pandasToDouble(const std::string& text, double& value) {
  const std::array<double, 309>& e = powersOfTen();
  const char* p = text.c_str();
  while (isSpace(*p)) {
    p++;
  }
  bool negative = false;
  if (*p == '-' || *p == '+') {
    negative = *p == '-';
    p++;
  }
  double number = 0.;
  int exponent = 0;
  int num_digits = 0;
  int num_decimals = 0;
  const int max_digits = 17;
  while (isDigit(*p)) {
    if (num_digits < max_digits) {
      number = number * 10. + (*p - '0');
      num_digits++;
    } else {
      ++exponent;
    }
    p++;
  }
  if (*p == '.') {
    p++;
    while (num_digits < max_digits && isDigit(*p)) {
      number = number * 10. + (*p - '0');
      p++;
      num_digits++;
      num_decimals++;
    }
    if (num_digits >= max_digits) {
      while (isDigit(*p)) {
        ++p;
      }
    }
    exponent -= num_decimals;
  }
  bool ok = num_digits > 0;
  if (ok) {
    if (negative) {
      number = -number;
    }
    if (*p == 'e' || *p == 'E') {
      bool exp_negative = false;
      ++p;
      if (*p == '-' || *p == '+') {
        exp_negative = *p == '-';
        p++;
      }
      int exp_digits = 0;
      int n = 0;
      while (exp_digits < max_digits && isDigit(*p)) {
        n = n * 10 + (*p - '0');
        exp_digits++;
        p++;
      }
      exponent = exp_negative ? exponent - n : exponent + n;
      // no digits after the 'e': un-consume it
      if (exp_digits == 0) {
        p--;
      }
    }
    if (exponent > 308) {
      ok = false;
    } else if (exponent > 0) {
      number *= e[exponent];
    } else if (exponent < -308) {
      if (exponent < -616) {
        number = 0.;
      } else {
        number /= e[-308 - exponent];
        number /= e[308];
      }
    } else {
      number /= e[-exponent];
    }
    if (std::isinf(number)) {
      ok = false;
    }
    while (isSpace(*p)) {
      p++;
    }
    ok = ok && *p == '\0';
  }
  if (ok) {
    value = number;
    return true;
  }
  std::string lower;
  for (char c : text) {
    lower += static_cast<char>(c >= 'A' && c <= 'Z' ? c - 'A' + 'a' : c);
  }
  if (lower == "inf" || lower == "+inf" || lower == "infinity" ||
      lower == "+infinity") {
    value = std::numeric_limits<double>::infinity();
    return true;
  }
  if (lower == "-inf" || lower == "-infinity") {
    value = -std::numeric_limits<double>::infinity();
    return true;
  }
  return false;
}

// ------------------------------------------------------------ PyRandom
PyRandom::PyRandom(uint64_t seed) {
  // random.seed(int): the absolute value as 32-bit words, least
  // significant first ([0] for zero)
  std::vector<uint32_t> key;
  while (seed > 0) {
    key.push_back(static_cast<uint32_t>(seed & 0xffffffffU));
    seed >>= 32;
  }
  if (key.empty()) {
    key.push_back(0);
  }
  initByArray(key);
}

void PyRandom::initGenrand(uint32_t s) {
  mt_[0] = s;
  for (int i = 1; i < kN; i++) {
    mt_[i] = 1812433253U * (mt_[i - 1] ^ (mt_[i - 1] >> 30)) +
             static_cast<uint32_t>(i);
  }
  index_ = kN;
}

void PyRandom::initByArray(const std::vector<uint32_t>& key) {
  initGenrand(19650218U);
  int i = 1;
  size_t j = 0;
  int key_length = static_cast<int>(key.size());
  for (int k = kN > key_length ? kN : key_length; k > 0; k--) {
    mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1664525U)) +
             key[j] + static_cast<uint32_t>(j);
    i++;
    j++;
    if (i >= kN) {
      mt_[0] = mt_[kN - 1];
      i = 1;
    }
    if (j >= key.size()) {
      j = 0;
    }
  }
  for (int k = kN - 1; k > 0; k--) {
    mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1566083941U)) -
             static_cast<uint32_t>(i);
    i++;
    if (i >= kN) {
      mt_[0] = mt_[kN - 1];
      i = 1;
    }
  }
  mt_[0] = 0x80000000U;
}

uint32_t PyRandom::genrandUint32() {
  static const uint32_t mag01[2] = {0x0U, 0x9908b0dfU};
  const uint32_t upper = 0x80000000U;
  const uint32_t lower = 0x7fffffffU;
  uint32_t y;
  if (index_ >= kN) {
    int kk = 0;
    for (; kk < kN - kM; kk++) {
      y = (mt_[kk] & upper) | (mt_[kk + 1] & lower);
      mt_[kk] = mt_[kk + kM] ^ (y >> 1) ^ mag01[y & 0x1U];
    }
    for (; kk < kN - 1; kk++) {
      y = (mt_[kk] & upper) | (mt_[kk + 1] & lower);
      mt_[kk] = mt_[kk + (kM - kN)] ^ (y >> 1) ^ mag01[y & 0x1U];
    }
    y = (mt_[kN - 1] & upper) | (mt_[0] & lower);
    mt_[kN - 1] = mt_[kM - 1] ^ (y >> 1) ^ mag01[y & 0x1U];
    index_ = 0;
  }
  y = mt_[index_++];
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680U;
  y ^= (y << 15) & 0xefc60000U;
  y ^= (y >> 18);
  return y;
}

int PyRandom::randbelow(int n) {
  int k = 0;
  while ((n >> k) > 0) {
    k++;
  }
  uint32_t r = getrandbits(k);
  while (r >= static_cast<uint32_t>(n)) {
    r = getrandbits(k);
  }
  return static_cast<int>(r);
}

// random.sample's list branch (n <= setsize = 21 for k <= 5)
std::vector<int> PyRandom::sample(int n, int k) {
  std::vector<int> pool(n);
  for (int i = 0; i < n; i++) {
    pool[i] = i;
  }
  std::vector<int> result(k);
  for (int i = 0; i < k; i++) {
    int j = randbelow(n - i);
    result[i] = pool[j];
    pool[j] = pool[n - i - 1];
  }
  return result;
}

}  // namespace lib_numeric

}  // namespace toplib

}  // namespace toppic
