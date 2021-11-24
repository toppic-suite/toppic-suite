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

#include <catch.hpp>

#include "common/base/support_peak_type_base.hpp"

namespace toppic {

TEST_CASE("support peak type initializing"){
    SPTypeBase::initBase();

    SECTION("checking data at support peak type address"){
        REQUIRE(SPTypeBase::getSPTypePtrByName("N_TERM")->getId() == 0);
        REQUIRE(SPTypeBase::getSPTypePtrById(1)->getName() == "C_TERM");
        REQUIRE(SPTypeBase::getSPTypePtrByName("N_TERM")->getXmlElementName() == "support_peak_type");
    }
    //invalid pointers
    REQUIRE(SPTypeBase::getSPTypePtrByName(" ") == nullptr);
    REQUIRE(SPTypeBase::getSPTypePtrById(100) == nullptr);

}

}