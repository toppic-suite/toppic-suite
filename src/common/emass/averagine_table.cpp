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

#include "common/emass/averagine_table.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <vector>

#include "common/emass/emass.hpp"
#include "common/util/logger.hpp"

namespace toppic {

namespace emass {

namespace {

// EnvBase::getDistrEntryNum() and getDistrMassInterval()
const int kEntryNum = 11000;
const double kMassInterval = 10.0;
// emass's print limit, relative to the apex
const double kTailLimit = 5e-9;

void writeEntry(std::ostream &out, double mass,
                const Composition &composition, int sulphur_num) {
  Formula formula = averagineFormulaWithSulphur(
      mass, sulphur_num, AtomRounding::kRoundEach, composition);
  IsotopeDistribution distr = getDistribution(formula);
  const std::vector<double> &rel = distr.rel_abundances;
  size_t apex = std::distance(rel.begin(),
                              std::max_element(rel.begin(), rel.end()));
  size_t end = apex + 1;
  while (end < rel.size() && rel[end] >= kTailLimit) {
    end++;
  }
  if (std::any_of(rel.begin(), rel.begin() + end,
                  [](double x) { return x <= 0; })) {
    LOG_ERROR("Non-positive isotope abundance at " << mass << " Da.");
    exit(EXIT_FAILURE);
  }
  char buf[64];
  out << "formula: ";
  for (int e = 0; e < kElementNum; e++) {
    if (formula[e] != 0) {
      out << getSymbol(static_cast<Element>(e))
          << static_cast<long>(formula[e]);
    }
  }
  std::snprintf(buf, sizeof(buf), "%.6f", distr.mono_mass);
  out << " charge: 0 limit: 0.000000 mass: " << buf << "\n";
  for (size_t k = 0; k < end; k++) {
    std::snprintf(buf, sizeof(buf), "%.6f %.7g",
                  distr.mono_mass + distr.mass_offsets[k], rel[k] * 100);
    out << buf << "\n";
  }
}

}  // namespace

void writeAveragineTable(std::ostream &out, const Composition &composition,
                         int sulphur_num) {
  for (int i = 0; i < kEntryNum; i++) {
    if (i > 0) {
      out << "\n";
    }
    writeEntry(out, kMassInterval * std::max(i, 1), composition,
               sulphur_num);
  }
}

void writeAveragineTable(const std::string &file_name,
                         const Composition &composition, int sulphur_num) {
  std::ofstream out(file_name);
  if (!out.is_open()) {
    LOG_ERROR("Cannot open the averagine table file " << file_name << "!");
    exit(EXIT_FAILURE);
  }
  writeAveragineTable(out, composition, sulphur_num);
  out.close();
  if (out.fail()) {
    LOG_ERROR("Cannot write the averagine table file " << file_name << "!");
    exit(EXIT_FAILURE);
  }
}

}  // namespace emass

}  // namespace toppic
