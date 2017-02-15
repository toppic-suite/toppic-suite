// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_MASS_CONSTANT_HPP_
#define PROT_MASS_CONSTANT_HPP_

namespace prot {

class MassConstant {
 public:
  /**
   * Returns the mass of an ammonia (NH3).
   */
  static double getAmmoniaMass() {return 17.0265; }

  /**
   * Returns the mass difference between two isotopic peaks.
   */
  static double getIsotopeMass() { /* from Thrash paper */ return 1.00235; }

  /**
   * Returns the mass of an oxygen molecular.
   */
  static double getOxygenMass() {return 15.9949; }

  /**
   * Returns the mass of a proton. 
   */
  static double getProtonMass() {return 1.007276; }

  /**
   * Returns the mass of a water molecule (H2O).
   */
  static double getWaterMass() {return 18.010565; }

  /**
   * Returns the shift between the mass of a neutral y ion 
   * and its corresponding suffix residue mass. 
   */
  static double getYIonShift() {return 18.010565; }
};

}
#endif
