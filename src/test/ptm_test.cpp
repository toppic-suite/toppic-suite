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

#include "common/base/ptm_base.hpp"

namespace toppic{

TEST_CASE("ptm initialization"){
    PtmBase::initBase();//initialize amino acid base

    SECTION("get ptm pointer by abbreviated name"){
        //expected to return valid pointer
        PtmPtr none = PtmBase::getPtmPtrByAbbrName("No PTM");
        PtmPtr triiod = PtmBase::getPtmPtrByAbbrName("Triiodothyronine");
        PtmPtr sulfide = PtmBase::getPtmPtrByAbbrName("Sulfide");

        //examine other information at that pointer address
        SECTION("checking unimodID"){
            REQUIRE(triiod->getUnimodId() == 397);
        }
        SECTION("checking mass comparison"){
            REQUIRE(none->cmpMassInc(none, triiod));//comparing mass. true if a<b
            REQUIRE(none->cmpMassInc(sulfide, triiod));
        }
        SECTION("checking if abbr. name is same"){
            REQUIRE(sulfide->isSame(sulfide));
        }
        SECTION("checking getter of element name"){
            REQUIRE(triiod->getXmlElementName() == "ptm");
        }

        SECTION ("checking getter of ptm pointer (getPtmPtr)"){
            REQUIRE(PtmBase::getPtmPtr(none) == none);
            REQUIRE(PtmBase::getPtmPtr(triiod) == triiod);
        }
    
        //expected to return null pointer
        REQUIRE(PtmBase::getPtmPtrByAbbrName(" ") == nullptr);
        REQUIRE(PtmBase::getPtmPtrByAbbrName("PREC") == nullptr);//this is ion name
    }
    SECTION("Checks if the list contains an amino acid with the specific name"){
        REQUIRE(PtmBase::containsAbbrName("TPQ"));
        REQUIRE(PtmBase::containsAbbrName("nesyl") == false);//part of Farnesyl
    }

}

}


