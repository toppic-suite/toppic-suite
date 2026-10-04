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

#ifndef TOPPIC_COMMON_EMASS_AVERAGINE_TABLE_HPP_
#define TOPPIC_COMMON_EMASS_AVERAGINE_TABLE_HPP_

#include <ostream>
#include <string>

#include "common/emass/isotope_data.hpp"

namespace toppic {

namespace emass {

// TopFD's averagine table res/base_data/theo_patt.txt, in the format EnvBase
// reads (src/ms/env/env_base.cpp, env.cpp): 11,000 entries, one per 10 Da of
// monoisotopic mass, each a formula line followed by "mass intensity" lines
// (intensity in % of the apex) and a blank line between entries. Entry i is
// the averagine of the monoisotopic mass 10 i Da with each of the five atom
// numbers rounded (AtomRounding::kRoundEach; entry 0 repeats entry 1). Every
// peak from the monoisotopic one is listed, since Env::Env takes the first
// listed peak as the monoisotopic one; the list ends after the apex at the
// last peak of at least 5e-9 of the apex (emass's print limit). Masses have
// six decimals, intensities seven significant digits. With the default
// human averagine the output is the table of toppic_test experiment 111.
void writeAveragineTable(std::ostream &out,
                         const Composition &composition = kHumanAveragine);

void writeAveragineTable(const std::string &file_name,
                         const Composition &composition = kHumanAveragine);

}  // namespace emass

}  // namespace toppic

#endif
