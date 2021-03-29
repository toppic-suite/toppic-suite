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

#include "common/base/prot_mod_base.hpp"

namespace toppic{

TEST_CASE("prot_mod initialization"){
    ProtModBase::initBase();//initialize mod base

    SECTION("get prot_mod pointer by name"){
        //expected to return valid pointer
        ProtModPtr nme_ace_ser = ProtModBase::getProtModPtrByName("NME_ACETYLATION_SER");
        ProtModPtr none = ProtModBase::getProtModPtrByName("NONE");

        REQUIRE(nme_ace_ser != ProtModPtr(nullptr));
        REQUIRE(none != ProtModPtr(nullptr));

        SECTION("get prot mod name and type"){
            REQUIRE(nme_ace_ser->getName() == "NME_ACETYLATION_SER");
            REQUIRE(none->getType() == "NONE");
        }

        SECTION("get pointers for mod and trunc"){
            REQUIRE(nme_ace_ser->getProtShift() != NULL);
            REQUIRE(nme_ace_ser->getPepShift() != NULL);
        }

        SECTION ("checking isAcetylation"){
            REQUIRE(nme_ace_ser->isAcetylation());
            REQUIRE(none->isAcetylation() == false);
        }
        SECTION("checking prot mod xml element name"){
            REQUIRE(nme_ace_ser->getXmlElementName() == "prot_mod");
        }
        //expected to return null pointer
        REQUIRE(ProtModBase::getProtModPtrByName(" ") == nullptr);
        REQUIRE(ProtModBase::getProtModPtrByName("Threon") == nullptr);//this is part of Threonine, which is valid name
    }

    SECTION("get prot_mod pointer by type"){
        //expected to return valid pointer
        REQUIRE(ProtModBase::getProtModPtrByType("NME").size() > 0 );//check the return pointer vector is not null
        REQUIRE(ProtModBase::getProtModPtrByType("M_ACETYLATION").size() > 0); 
        //expected to return null pointer
        REQUIRE(ProtModBase::getProtModPtrByType("ABC").size() < 1); 
        REQUIRE(ProtModBase::getProtModPtrByType(" ").size() < 1); 
    }


    }

}

