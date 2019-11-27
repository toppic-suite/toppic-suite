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