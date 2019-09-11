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

#include "common/base/residue.hpp"
#include "common/base/residue_freq.hpp"

namespace toppic {

ResidueFreq::ResidueFreq(AminoAcidPtr acid_ptr, 
                         PtmPtr ptm_ptr, double freq):
    Residue(acid_ptr, ptm_ptr),
    freq_(freq) {}

}  // namespace toppic
