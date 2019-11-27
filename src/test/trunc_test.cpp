
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