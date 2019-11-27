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