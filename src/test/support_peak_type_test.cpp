
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