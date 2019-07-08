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

#ifndef TOPPIC_COMMON_BASE_ION_TYPE_DATA_HPP_
#define TOPPIC_COMMON_BASE_ION_TYPE_DATA_HPP_

# include <string>

namespace toppic {

std::string ion_type_base_data = R"(
<ion_type_list>
<ion_type>
	<name>B</name>
	<n_term>true</n_term>
	<shift>0</shift>
</ion_type>
<ion_type>
	<name>Y</name>
	<n_term>false</n_term>
	<shift>18.0106</shift>
</ion_type>
<ion_type>
	<name>C</name>
	<n_term>true</n_term>
	<shift>17.0265</shift>
</ion_type>
<ion_type>
	<name>Z_DOT</name>
	<n_term>false</n_term>
	<shift>1.9919</shift>
</ion_type>
<ion_type>
	<name>A</name>
	<n_term>true</n_term>
	<shift>-26.9871</shift>
</ion_type>
<ion_type>
	<name>X</name>
	<n_term>false</n_term>
	<shift>44.9977</shift>
</ion_type>
<ion_type>
	<name>PREC</name>
	<n_term>false</n_term>
	<shift>18.0106</shift>
</ion_type>
</ion_type_list>
)";

} 

#endif
