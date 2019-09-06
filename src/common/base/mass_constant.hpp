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


#ifndef TOPPIC_COMMON_BASE_MASS_CONSTANT_HPP_
#define TOPPIC_COMMON_BASE_MASS_CONSTANT_HPP_

namespace toppic {

// mass values are obtained from protpi.ch 
namespace mass_constant {
/**
 * Returns the mass of an ammonia (NH3).
 */
inline double getAmmoniaMass() {return 17.026549101073; }

/**
 * Returns the mass difference between two isotopic peaks.
 * The value is from the paper of the deconvolution tool Thrash.
 * We also double checked with the averagine molecule C494 H776 N136 O148 S4
 */
inline double getIsotopeMass() { /* from Thrash paper */ return 1.00235; }

/**
 * Returns the mass of an oxygen molecular.
 */
inline double getOxygenMass() {return 15.9949146195616; }

/**
 * Returns the mass of a proton. 
 */
inline double getProtonMass() {/* from wikipedia */ return 1.007276466879; }

/**
 * Returns the mass of a water molecule (H2O).
 */
inline double getWaterMass() {return 18.0105646837036; }

/**
 * Returns the shift between the mass of a neutral y ion 
 * and its corresponding suffix residue mass. 
 */
inline double getYIonShift() {return 18.0105646837036; }

}  // namespace mass_constant

}  // namespace toppic

#endif
