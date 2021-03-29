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
