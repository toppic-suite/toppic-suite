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
