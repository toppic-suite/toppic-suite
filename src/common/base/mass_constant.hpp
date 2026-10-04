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

#ifndef TOPPIC_COMMON_BASE_MASS_CONSTANT_HPP_
#define TOPPIC_COMMON_BASE_MASS_CONSTANT_HPP_

namespace toppic {

// Monoisotopic masses are sums of the AME2020 atomic masses (the Atomic Mass
// Evaluation 2020, the masses IUPAC's CIAAW uses): 1H 1.007825031898,
// 14N 14.00307400425, 16O 15.99491461926. The proton mass is CODATA 2022's.
namespace mass_constant {
/**
 * Returns the mass of an ammonia (NH3).
 */
inline double getAmmoniaMass() { /* AME2020 */ return 17.026549099944; }

/**
 * Returns the mass difference between two isotopic peaks.
 * The value is from the paper of the deconvolution tool Thrash.
 * We also double checked with the averagine molecule C494 H776 N136 O148 S4
 */
inline double getIsotopeMass() { /* from Thrash paper */ return 1.00235; }

/**
 * Returns the mass of an oxygen atom.
 */
inline double getOxygenMass() { /* AME2020 */ return 15.99491461926; }

/**
 * Returns the mass of a proton.
 */
inline double getProtonMass() { /* CODATA 2022 */ return 1.0072764665789; }

/**
 * Returns the mass of a water molecule (H2O).
 */
inline double getWaterMass() { /* AME2020 */ return 18.010564683056; }

/**
 * Returns the shift between the mass of a neutral y ion
 * and its corresponding suffix residue mass.
 */
inline double getYIonShift() { /* AME2020 */ return 18.010564683056; }

}  // namespace mass_constant

}  // namespace toppic

#endif
