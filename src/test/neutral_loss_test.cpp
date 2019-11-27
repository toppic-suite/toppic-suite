#include <string>
#include <catch.hpp>

#include "common/base/neutral_loss_base.hpp"

using namespace toppic;

namespace toppic{

TEST_CASE("neutral loss initialization"){
    NeutralLossBase::initBase();//initialize ion base

    SECTION("get neutral loss pointer by name"){
        //expected to return valid pointer
        NeutralLossPtr water = NeutralLossBase::getNeutralLossPtrByName("Water");
        NeutralLossPtr none = NeutralLossBase::getNeutralLossPtrByName("NONE");

        //examine other information such as composition at that pointer address
        SECTION("examine neutral loss pointer address"){
            REQUIRE(none->getMass() == 0.0);
            REQUIRE(water->getName() == "Water");
            REQUIRE(water->getXmlElementName() == "neutral_loss");
        
        }

        //expected to return null pointer
        REQUIRE(NeutralLossBase::getNeutralLossPtrByName(" ") == nullptr);
        REQUIRE(NeutralLossBase::getNeutralLossPtrByName("Glutamine") == nullptr);//this is amino acid name
    }

}

}
