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

#ifndef TOPPIC_COMMON_BASE_NEUTRAL_LOSS_DATA_HPP_
#define TOPPIC_COMMON_BASE_NEUTRAL_LOSS_DATA_HPP_

# include <string>

namespace toppic {

std::string neutral_loss_base_data = R"(
<neutral_loss_list>
<neutral_loss>
	<name>NONE</name>
	<mass>0</mass>
</neutral_loss>
<neutral_loss>
	<name>Water</name>
	<mass>18.0106</mass>
</neutral_loss>
<neutral_loss>
	<name>Ammonia</name>
	<mass>17.0265</mass>
</neutral_loss>
</neutral_loss_list>
)";

} 

#endif
