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

#ifndef TOPPIC_COMMON_BASE_PROT_MOD_DATA_HPP_
#define TOPPIC_COMMON_BASE_PROT_MOD_DATA_HPP_

# include <string>

namespace toppic {

std::string prot_mod_base_data = R"(
<prot_mod_list>
  <prot_mod>
    <name>NONE</name>
    <type>NONE</type>
    <truncation>
      <name>NONE</name>
    </truncation>
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
  </prot_mod>
  <prot_mod>
    <name>M_ACETYLATION</name>
    <type>M_ACETYLATION</type>
    <truncation>
      <name>NONE</name>
    </truncation>
    <mod>
      <ori_residue>
        <amino_acid>
	        <name>Methionine</name>
        </amino_acid>
        <ptm>
          <abbreviation>No PTM</abbreviation>
        </ptm>
      </ori_residue>
      <mod_residue>
        <amino_acid>
	        <name>Methionine</name>
        </amino_acid>
        <ptm>
          <abbreviation>Acetyl</abbreviation>
        </ptm>
      </mod_residue>
    </mod>
  </prot_mod>
  <prot_mod>
    <name>NME</name>
    <type>NME</type>
    <truncation>
      <name>NME</name>
    </truncation>
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
  </prot_mod>
  <prot_mod>
    <name>NME_ACETYLATION_ALA</name>
    <type>NME_ACETYLATION</type>
    <truncation>
      <name>NME</name>
    </truncation>
    <mod>
      <ori_residue>
        <amino_acid>
	        <name>Alanine</name>
        </amino_acid>
        <ptm>
          <abbreviation>No PTM</abbreviation>
        </ptm>
      </ori_residue>
      <mod_residue>
        <amino_acid>
	        <name>Alanine</name>
        </amino_acid>
        <ptm>
          <abbreviation>Acetyl</abbreviation>
        </ptm>
      </mod_residue>
    </mod>
  </prot_mod>
  <prot_mod>
    <name>NME_ACETYLATION_GLY</name>
    <type>NME_ACETYLATION</type>
    <truncation>
      <name>NME</name>
    </truncation>
    <mod>
      <ori_residue>
        <amino_acid>
	        <name>Glycine</name>
        </amino_acid>
        <ptm>
          <abbreviation>No PTM</abbreviation>
        </ptm>
      </ori_residue>
      <mod_residue>
        <amino_acid>
	        <name>Glycine</name>
        </amino_acid>
        <ptm>
          <abbreviation>Acetyl</abbreviation>
        </ptm>
      </mod_residue>
    </mod>
  </prot_mod>
  <prot_mod>
    <name>NME_ACETYLATION_SER</name>
    <type>NME_ACETYLATION</type>
    <truncation>
      <name>NME</name>
    </truncation>
    <mod>
      <ori_residue>
        <amino_acid>
	        <name>Serine</name>
        </amino_acid>
        <ptm>
          <abbreviation>No PTM</abbreviation>
        </ptm>
      </ori_residue>
      <mod_residue>
        <amino_acid>
	        <name>Serine</name>
        </amino_acid>
        <ptm>
          <abbreviation>Acetyl</abbreviation>
        </ptm>
      </mod_residue>
    </mod>
  </prot_mod>
  <prot_mod>
    <name>NME_ACETYLATION_THR</name>
    <type>NME_ACETYLATION</type>
    <truncation>
      <name>NME</name>
    </truncation>
    <mod>
      <ori_residue>
        <amino_acid>
	        <name>Threonine</name>
        </amino_acid>
        <ptm>
          <abbreviation>No PTM</abbreviation>
        </ptm>
      </ori_residue>
      <mod_residue>
        <amino_acid>
	        <name>Threonine</name>
        </amino_acid>
        <ptm>
          <abbreviation>Acetyl</abbreviation>
        </ptm>
      </mod_residue>
    </mod>
  </prot_mod>
  <prot_mod>
    <name>NME_ACETYLATION_VAL</name>
    <type>NME_ACETYLATION</type>
    <truncation>
      <name>NME</name>
    </truncation>
    <mod>
      <ori_residue>
        <amino_acid>
	        <name>Valine</name>
        </amino_acid>
        <ptm>
          <abbreviation>No PTM</abbreviation>
        </ptm>
      </ori_residue>
      <mod_residue>
        <amino_acid>
	        <name>Valine</name>
        </amino_acid>
        <ptm>
          <abbreviation>Acetyl</abbreviation>
        </ptm>
      </mod_residue>
    </mod>
  </prot_mod>
  <prot_mod>
    <name>NME_ACETYLATION_CYS</name>
    <type>NME_ACETYLATION</type>
    <truncation>
      <name>NME</name>
    </truncation>
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
          <abbreviation>Acetyl</abbreviation>
        </ptm>
      </mod_residue>
    </mod>
  </prot_mod>
</prot_mod_list>
)";

} 

#endif
