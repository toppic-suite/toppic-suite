#include <catch.hpp>

#include "common/base/mod_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/amino_acid_base.hpp"
#include "common/base/mod_util.hpp"

namespace toppic{
namespace mod_util {
    TEST_CASE("mod initializing"){
        ModBase::initBase();
        AminoAcidPtr cysteine = AminoAcidBase::getAminoAcidPtrByName("Cysteine");
        AminoAcidPtr isol = AminoAcidBase::getAminoAcidPtrByName("Isoleucine");

        ResiduePtr residuePtr_C = ResidueBase::getBaseResiduePtr(cysteine);
        ResiduePtr residuePtr_I = ResidueBase::getBaseResiduePtr(isol);

        ModPtr modPtr =  ModBase::getBaseModPtr(residuePtr_C, residuePtr_C);
        ResiduePtrVec residuePtrVec;
        
        residuePtrVec.push_back(residuePtr_I);

        REQUIRE(modPtr->getOriResiduePtr() == residuePtr_C);
        REQUIRE(modPtr->getModResiduePtr() == residuePtr_C);
        REQUIRE(modPtr->isSame(modPtr));
        REQUIRE(modPtr->getXmlElementName() == "mod");
        REQUIRE(modPtr->getShift() == 0);

        SECTION("testing mod utils"){
            REQUIRE(geneResidueListWithMod(residuePtrVec, geneFixedModList("C57")).size() > 0);
            REQUIRE(getModMassVec(geneFixedModList("C58")).size() > 0);
        }
    }
}
}