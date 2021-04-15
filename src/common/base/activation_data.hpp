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

#ifndef TOPPIC_COMMON_BASE_ACTIVATION_DATA_HPP_
#define TOPPIC_COMMON_BASE_ACTIVATION_DATA_HPP_

# include <string>

namespace toppic {

std::string activation_base_data = R"(
<activation_list>
<activation>
	<name>CID</name>
	<n_ion_type>B</n_ion_type>
	<c_ion_type>Y</c_ion_type>
</activation>
<activation>
	<name>HCD</name>
	<n_ion_type>B</n_ion_type>
	<c_ion_type>Y</c_ion_type>
</activation>
<activation>
	<name>ETD</name>
	<n_ion_type>C</n_ion_type>
	<c_ion_type>Z_DOT</c_ion_type>
</activation>
<activation>
	<name>MPD</name>
	<n_ion_type>A</n_ion_type>
	<c_ion_type>Y</c_ion_type>
</activation>
<activation>
	<name>UVPD</name>
	<n_ion_type>A</n_ion_type>
	<c_ion_type>Y</c_ion_type>
</activation>
</activation_list>
)";

} 

#endif
