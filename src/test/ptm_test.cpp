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


