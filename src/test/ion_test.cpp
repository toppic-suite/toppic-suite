#include <string>
#include <catch.hpp>

#include "common/base/ion_type_base.hpp"

namespace toppic{

TEST_CASE("ion initialization"){
    IonTypeBase::initBase();//initialize ion base

    SECTION("get ion type pointer by name"){
        //expected to return valid pointer
        IonTypePtr typeB = IonTypeBase::getIonTypePtrByName("B");
        IonTypePtr typeA = IonTypeBase::getIonTypePtrByName("A");

        //examine other information such as composition at that pointer address
        SECTION("examine ion type pointer address"){
            REQUIRE(typeB->getName() == "B");
            REQUIRE(typeB->getShift() == 0);

            REQUIRE(typeA->isNTerm());
            REQUIRE(typeA->getXmlElementName() == "ion_type");

            //charge and pos were also working correctly when tested with random value for getter functions
        }

        //expected to return null pointer
        REQUIRE(IonTypeBase::getIonTypePtrByName(" ") == nullptr);
        REQUIRE(IonTypeBase::getIonTypePtrByName("V") == nullptr);
        REQUIRE(IonTypeBase::getIonTypePtrByName("Cysteine") == nullptr);//this is amino acid name
    }

}

}
