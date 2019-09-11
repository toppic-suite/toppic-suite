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

#ifndef TOPPIC_PRSM_PRSM_ALGO_HPP_
#define TOPPIC_PRSM_PRSM_ALGO_HPP_

#include <vector>

namespace toppic {

namespace prsm_algo {

// used in find matched mass pairs
bool increaseIJ(size_t i, size_t j, double deviation, double tolerance,
                const std::vector<double> &ms_masses,
                const std::vector<double> &theo_masses);

// compute ppos for ms_masses
std::vector<double> compMsMassPpos(const std::vector<double> &ms_masses,
                                   const std::vector<double> &theo_masses,
                                   double ppo);

// compute the number of matched theoretical masses (fragment ions)
double compNumMatchedTheoMasses(const std::vector<double> &ms_masses,
                                const std::vector<double> &theo_masses,
                                double ppo);

// compute the position of the last residue of a proteoform based
// on its n term shift
int getFirstResPos(double n_term_shift, const std::vector<double> &prm_masses);

// compute the position of the last residue of a proteoform based
// on its c term shift
int getLastResPos(double c_term_shift, const std::vector<double> &prm_masses);

}  // namespace prsm_algo

}  // namespace toppic

#endif
