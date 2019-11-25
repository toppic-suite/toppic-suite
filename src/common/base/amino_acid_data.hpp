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

#ifndef TOPPIC_COMMON_BASE_AMINO_ACID_DATA_HPP_
#define TOPPIC_COMMON_BASE_AMINO_ACID_DATA_HPP_

# include <string>

namespace toppic {

// Monoisotopic masses are computed using the IUPAC 
// atomic weight (version 2017). Average masses are 
// obtained from 
// https://proteomicsresource.washington.edu/protocols06/masses.php 

std::string amino_acid_base_data = R"(
<amino_acid_list>
<amino_acid>
	<name>None</name>
	<one_letter>X</one_letter>
	<three_letter>Xxx</three_letter>
	<composition></composition>
	<mono_mass>0.0</mono_mass>
	<average_mass>0.0</average_mass>
	<frequency>0.0</frequency>
</amino_acid>
<amino_acid>
	<name>Alanine</name>
	<one_letter>A</one_letter>
	<three_letter>Ala</three_letter>
	<composition>C3H5ON</composition>
	<mono_mass>71.0371137846</mono_mass>
	<average_mass>71.0779</average_mass>
	<frequency>0.0828</frequency>
</amino_acid>
<amino_acid>
	<name>Arginine</name>
	<one_letter>R</one_letter>
	<three_letter>Arg</three_letter>
	<composition>C6H12ON4</composition>
	<mono_mass>156.10111102302</mono_mass>
	<average_mass>156.18568</average_mass>
	<frequency>0.0554</frequency>
</amino_acid>
<amino_acid>
	<name>Asparagine</name>
	<one_letter>N</one_letter>
	<three_letter>Asn</three_letter>
	<composition>C4H6O2N2</composition>
	<mono_mass>114.04292744016</mono_mass>
	<average_mass>114.10264</average_mass>
	<frequency>0.0404</frequency>
</amino_acid>
<amino_acid>
	<name>Aspartic_acid</name>
	<one_letter>D</one_letter>
	<three_letter>Asp</three_letter>
	<composition>C4H5O3N</composition>
	<mono_mass>115.02694302280</mono_mass>
	<average_mass>115.0874</average_mass>
	<frequency>0.0545</frequency>
</amino_acid>
<amino_acid>
	<name>Cysteine</name>
	<one_letter>C</one_letter>
	<three_letter>Cys</three_letter>
	<composition>C3H5ONS</composition>
	<mono_mass>103.0091849595</mono_mass>
	<average_mass>103.1429</average_mass>
	<frequency>0.0136</frequency>
</amino_acid>
<amino_acid>
	<name>Glutamic_acid</name>
	<one_letter>E</one_letter>
	<three_letter>Glu</three_letter>
	<composition>C5H7O3N</composition>
	<mono_mass>129.04259308732</mono_mass>
	<average_mass>129.11398</average_mass>
	<frequency>0.0394</frequency>
</amino_acid>
<amino_acid>
	<name>Glutamine</name>
	<one_letter>Q</one_letter>
	<three_letter>Gln</three_letter>
	<composition>C5H8O2N2</composition>
	<mono_mass>128.05857750468</mono_mass>
	<average_mass>128.12922</average_mass>
	<frequency>0.0677</frequency>
</amino_acid>
<amino_acid>
	<name>Glycine</name>
	<one_letter>G</one_letter>
	<three_letter>Gly</three_letter>
	<composition>C2H3ON</composition>
	<mono_mass>57.02146372008</mono_mass>
	<average_mass>57.05132</average_mass>
	<frequency>0.0709</frequency>
</amino_acid>
<amino_acid>
	<name>Histidine</name>
	<one_letter>H</one_letter>
	<three_letter>His</three_letter>
	<composition>C6H7ON3</composition>
	<mono_mass>137.05891185752</mono_mass>
	<average_mass>137.13928</average_mass>
	<frequency>0.0227</frequency>
</amino_acid>
<amino_acid>
	<name>Isoleucine</name>
	<one_letter>I</one_letter>
	<three_letter>Ile</three_letter>
	<composition>C6H11ON</composition>
	<mono_mass>113.08406397816</mono_mass>
	<average_mass>113.15764</average_mass>
	<frequency>0.0599</frequency>
</amino_acid>
<amino_acid>
	<name>Leucine</name>
	<one_letter>L</one_letter>
	<three_letter>Leu</three_letter>
	<composition>C6H11ON</composition>
	<mono_mass>113.08406397816</mono_mass>
	<average_mass>113.15764</average_mass>
	<frequency>0.0967</frequency>
</amino_acid>
<amino_acid>
	<name>Lysine</name>
	<one_letter>K</one_letter>
	<three_letter>Lys</three_letter>
	<composition>C6H12ON2</composition>
	<mono_mass>128.09496301462</mono_mass>
	<average_mass>128.17228</average_mass>
	<frequency>0.0586</frequency>
</amino_acid>
<amino_acid>
	<name>Methionine</name>
	<one_letter>M</one_letter>
	<three_letter>Met</three_letter>
	<composition>C5H9ONS</composition>
	<mono_mass>131.04048508854</mono_mass>
	<average_mass>131.19606</average_mass>
	<frequency>0.0243</frequency>
</amino_acid>
<amino_acid>
	<name>Phenylalanine</name>
	<one_letter>F</one_letter>
	<three_letter>Phe</three_letter>
	<composition>C9H9ON</composition>
	<mono_mass>147.06841391364</mono_mass>
	<average_mass>147.17386</average_mass>
	<frequency>0.0386</frequency>
</amino_acid>
<amino_acid>
	<name>Proline</name>
	<one_letter>P</one_letter>
	<three_letter>Pro</three_letter>
	<composition>C5H7ON</composition>
	<mono_mass>97.05276384912</mono_mass>
	<average_mass>97.11518</average_mass>
	<frequency>0.0468</frequency>
</amino_acid>
<amino_acid>
	<name>Serine</name>
	<one_letter>S</one_letter>
	<three_letter>Ser</three_letter>
	<composition>C3H5O2N</composition>
	<mono_mass>87.03202840370</mono_mass>
	<average_mass>87.0773</average_mass>
	<frequency>0.0649</frequency>
</amino_acid>
<amino_acid>
	<name>Threonine</name>
	<one_letter>T</one_letter>
	<three_letter>Thr</three_letter>
	<composition>C4H7O2N</composition>
	<mono_mass>101.04767846822</mono_mass>
	<average_mass>101.10388</average_mass>
	<frequency>0.0532</frequency>
</amino_acid>
<amino_acid>
	<name>Tryptophan</name>
	<one_letter>W</one_letter>
	<three_letter>Trp</three_letter>
	<composition>C11H10ON2</composition>
	<mono_mass>186.07931295010</mono_mass>
	<average_mass>186.2099</average_mass>
	<frequency>0.0108</frequency>
</amino_acid>
<amino_acid>
	<name>Tyrosine</name>
	<one_letter>Y</one_letter>
	<three_letter>Tyr</three_letter>
	<composition>C9H9NO2</composition>
	<mono_mass>163.06332853274</mono_mass>
	<average_mass>163.17326</average_mass>
	<frequency>0.0291</frequency>
</amino_acid>
<amino_acid>
	<name>Valine</name>
	<one_letter>V</one_letter>
	<three_letter>Val</three_letter>
	<composition>C5H9NO</composition>
	<mono_mass>99.06841391364</mono_mass>
	<average_mass>99.13106</average_mass>
	<frequency>0.0688</frequency>
</amino_acid>
<amino_acid>
	<name>Selenocysteine</name>
	<one_letter>U</one_letter>
	<three_letter>Sec</three_letter>
    <composition>C3H5NOSe</composition>
	<mono_mass>150.9536363846</mono_mass>
	<average_mass>150.0379</average_mass>
	<frequency>0</frequency>
</amino_acid>
<amino_acid>
	<name>Ornithine</name>
	<one_letter>O</one_letter>
	<three_letter>Orn</three_letter>
    <composition>C5H12N2O2</composition>
	<mono_mass>132.08987763372</mono_mass>
	<average_mass>132.16098</average_mass>
	<frequency>0</frequency>
</amino_acid>
</amino_acid_list>
)";

} 

#endif
