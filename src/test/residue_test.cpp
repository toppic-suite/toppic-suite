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

#include "common/base/residue_base.hpp"
#include "common/base/amino_acid_base.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_util.hpp"

namespace toppic{
namespace residue_util {
TEST_CASE("residue initialization"){
    ResidueBase::initBase();

    //to get residue pointer, needs ptm and amino pointer first
    AminoAcidPtr lysine = AminoAcidBase::getAminoAcidPtrByName("Lysine");
    PtmPtr sulfide = PtmBase::getPtmPtrByAbbrName("Sulfide");
    ResiduePtr residuePtr = ResidueBase::getBaseResiduePtr(lysine, sulfide);

    SECTION("getting residue pointer"){
        REQUIRE(ResidueBase::getBaseResiduePtr(lysine) != nullptr);
        
            SECTION("checking residue class at residue pointer address"){
                REQUIRE(residuePtr->getAminoAcidPtr() == lysine);
                REQUIRE(residuePtr->getPtmPtr() == sulfide);
                REQUIRE(residuePtr->getMass() == lysine->getMonoMass() + sulfide->getMonoMass());
                REQUIRE(residuePtr->isSame(residuePtr));
                REQUIRE(residuePtr->getXmlElementName() == "residue");
        }
    }       
    SECTION("testing residue_util.cpp"){
        
        REQUIRE(isValidResidue('H'));
        REQUIRE(isValidResidue(8) == false);
        REQUIRE(isValidResidue('*') == false);
        REQUIRE(isValidResidue(' ') == false);
        
        REQUIRE(replaceResidueLetter('B') == 'D');
        
        REQUIRE(convertStrToResiduePtrVec("ABCDE").size() == 5);
        REQUIRE(convertStrToResiduePtrVec("0H S3-").size() == 2); // it is converted to HS
        
        REQUIRE(compResiduePtrVecMass("H") == 137.058911858638);

        REQUIRE(findResidue(convertStrToResiduePtrVec("AKD "), residuePtr));
    }
}


}
}