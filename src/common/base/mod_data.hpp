//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_COMMON_BASE_MOD_DATA_HPP_
#define TOPPIC_COMMON_BASE_MOD_DATA_HPP_

#include <string>

namespace toppic {

std::string mod_base_data = R"(
<mod_list>
  <mod>
    <ori_residue>
      <amino_acid>
        <name>None</name>
      </amino_acid>
      <ptm>
        <abbreviation>No PTM</abbreviation>
      </ptm>
    </ori_residue>
    <mod_residue>
      <amino_acid>
        <name>None</name>
      </amino_acid>
      <ptm>
        <abbreviation>No PTM</abbreviation>
      </ptm>
    </mod_residue>
  </mod>
  <mod>
    <ori_residue>
      <amino_acid>
	      <name>Cysteine</name>
      </amino_acid>
      <ptm>
        <abbreviation>No PTM</abbreviation>
      </ptm>
    </ori_residue>
    <mod_residue>
      <amino_acid>
	      <name>Cysteine</name>
      </amino_acid>
      <ptm>
        <abbreviation>Carbamidomethylation</abbreviation>
      </ptm>
    </mod_residue>
  </mod>
  <mod>
    <ori_residue>
      <amino_acid>
	      <name>Cysteine</name>
      </amino_acid>
      <ptm>
        <abbreviation>No PTM</abbreviation>
      </ptm>
    </ori_residue>
    <mod_residue>
      <amino_acid>
	      <name>Cysteine</name>
      </amino_acid>
      <ptm>
        <abbreviation>Carboxymethyl</abbreviation>
      </ptm>
    </mod_residue>
  </mod>
</mod_list>
)";

} 

#endif
