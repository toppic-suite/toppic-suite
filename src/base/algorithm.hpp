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


#ifndef PROT_BASE_ALGORITHM_HPP_
#define PROT_BASE_ALGORITHM_HPP_

#include <cstddef>
#include <vector>

namespace prot {

// used in find matched mass pairs 
bool increaseIJ(size_t i, size_t j, double deviation, double tolerance, 
                const std::vector<double> &ms_masses, 
                const std::vector<double> &theo_masses);

// compute ppos for ms_masses 
std::vector<double> compMsMassPpos(const std::vector<double> &ms_masses, 
                                   const std::vector<double> &theo_masses, 
                                   double ppo);

// compute the number of matched theoretical masses (fragment ions) 
double compNumMatchedTheoMasses (const std::vector<double> &ms_masses, 
                                 const std::vector<double> &theo_masses, 
                                 double ppo); 

// compute the position of the last residue of a proteoform based 
// on its n term shift
int getFirstResPos(double n_term_shift, const std::vector<double> &prm_masses);


// compute the position of the last residue of a proteoform based 
// on its c term shift 
int getLastResPos(double c_term_shift, const std::vector<double> &prm_masses);

}

#endif
