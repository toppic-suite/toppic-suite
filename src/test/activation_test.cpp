#include <catch.hpp>

#include "common/base/activation_base.hpp"
#include "common/base/ion_type_base.hpp"
namespace toppic{

TEST_CASE("checking activation"){
    IonTypeBase::initBase();
    ActivationBase::initBase();

    SECTION("getting activation pointer by name"){
        REQUIRE(ActivationBase::getActivationPtrByName(" ") == ActivationPtr(nullptr));

        SECTION("checking data at activation pointer address"){
            REQUIRE(ActivationBase::getActivationPtrByName("UVPD")->getNIonTypePtr() != nullptr);  
            REQUIRE(ActivationBase::getActivationPtrByName("ETD")->getCIonTypePtr() != nullptr); 

            ////////////////////////////////////////////////////////////////////////////////////////////
            //have to look to find out how exactly b shift is calculated
            //REQUIRE(ActivationBase::getActivationPtrByName("ETD")->getCShift()== 1.9919);  
            //REQUIRE(ActivationBase::getActivationPtrByName("UVPD")->getNShift() == -26.9871);  
            /////////////////////////////////////////////////////////////////////////////////////////

            REQUIRE(ActivationBase::getActivationPtrByName("CID")->getName() == "CID");  
            REQUIRE(ActivationBase::getActivationPtrByName("HCD")->getXmlElementName() == "activation");
        }
    }
}

}
