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

#include "common/util/str_util.hpp"

namespace toppic {
namespace str_util {
    TEST_CASE("test string processing"){

        REQUIRE(split("amino, ptm, ace", ",").size() == 3);

        REQUIRE(toString(true) == "1");
        REQUIRE(toString(365) == "365");
        REQUIRE(toString(-0.12345) == "-1.2345000000e-01");
        REQUIRE(toString(-0.12345, 1) == "-0.1");
        REQUIRE(toScientificStr(-0.12345, 3) == "-1.23e-01");
        REQUIRE(scientificToDouble("1.314e+1") == 13.14);

    }
}
}