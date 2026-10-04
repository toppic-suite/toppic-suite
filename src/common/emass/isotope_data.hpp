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

#ifndef TOPPIC_COMMON_EMASS_ISOTOPE_DATA_HPP_
#define TOPPIC_COMMON_EMASS_ISOTOPE_DATA_HPP_

#include <array>
#include <string>
#include <vector>

namespace toppic {

namespace emass {

// The five elements of protein averagine, in the order of the formulas.
enum Element { kCarbon = 0, kHydrogen, kNitrogen, kOxygen, kSulfur,
  kElementNum };

// An isotope: its AME2020 atomic mass (the Atomic Mass Evaluation 2020, Wang
// et al., Chinese Phys. C 45:030003, 2021, as CIAAW lists it) and its IUPAC
// representative abundance (Rosman and Taylor 1998, as NIST lists them). An
// abundance of 0 fills a mass number without a stable isotope (35S).
struct Isotope {
  double mass;
  double abundance;
};

// The isotopes of an element, indexed by the number of extra neutrons from
// the lightest one.
const std::vector<Isotope>& getIsotopes(Element element);

// The element's symbol: C, H, N, O or S.
const std::string& getSymbol(Element element);

// Atoms per amino acid residue.
using Composition = std::array<double, kElementNum>;

// The human averagine: the mean residue of the 20 standard amino acids in
// the UniProt human reference proteome (UP000005640, release 2026_03),
// weighted by their frequencies (toppic_test experiment 108).
inline constexpr Composition kHumanAveragine = {4.920188, 7.764730, 1.374482,
                                                1.484151, 0.044276};

// The human averagine without sulphur: the mean residue of the same
// proteome over the 20 standard amino acids minus cysteine and methionine
// (the only ones with sulphur), weighted by their frequencies renormalized
// over the 18 (toppic_test experiment 112). C, H, N and O only.
inline constexpr Composition kHumanNoSulphurAveragine = {
    4.964606, 7.803735, 1.391831, 1.506580, 0.0};

// Senko's averagine (Senko, Beu and McLafferty, JASMS 6:229, 1995).
inline constexpr Composition kSenkoAveragine = {4.9384, 7.7583, 1.3577,
                                                1.4773, 0.0417};

}  // namespace emass

}  // namespace toppic

#endif
