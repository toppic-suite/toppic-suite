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

#include "common/base/trunc_base.hpp"

namespace toppic {

TEST_CASE("trunc initializing"){
    TruncBase::initBase();

    TruncPtr none = TruncBase::getTruncPtrByName("NONE");
    TruncPtr nme = TruncBase::getTruncPtrByName("NME");

    SECTION("checking data at trunc pointer address"){
        REQUIRE(nme->getName() == "NME");
        REQUIRE(nme->getTruncLen() == 1);
        REQUIRE(none->getXmlElementName() == "truncation");
    }
    //invalid pointers
    REQUIRE(TruncBase::getTruncPtrByName(" ") == nullptr);
    REQUIRE(TruncBase::getTruncPtrByName("Alanine") == nullptr);

}

}