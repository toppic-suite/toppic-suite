//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_BASE_MASS_CONSTANT_HPP_
#define TOPPIC_BASE_MASS_CONSTANT_HPP_

namespace toppic {

namespace mass_constant {
/**
 * Returns the mass of an ammonia (NH3).
 */
inline double getAmmoniaMass() {return 17.0265; }

/**
 * Returns the mass difference between two isotopic peaks.
 */
inline double getIsotopeMass() { /* from Thrash paper */ return 1.00235; }

/**
 * Returns the mass of an oxygen molecular.
 */
inline double getOxygenMass() {return 15.9949; }

/**
 * Returns the mass of a proton. 
 */
inline double getProtonMass() {return 1.007276; }

/**
 * Returns the mass of a water molecule (H2O).
 */
inline double getWaterMass() {return 18.010565; }

/**
 * Returns the shift between the mass of a neutral y ion 
 * and its corresponding suffix residue mass. 
 */
inline double getYIonShift() {return 18.010565; }

}  // namespace mass_constant

}  // namespace toppic

#endif
