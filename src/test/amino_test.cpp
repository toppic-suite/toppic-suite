#include <string>
#include <catch.hpp>

#include "common/base/amino_acid_base.hpp"

namespace toppic{

TEST_CASE("amino acid initialization"){
    AminoAcidBase::initBase();//initialize amino acid base

    SECTION("get amino acid pointer by name"){
        //expected to return valid pointer
        REQUIRE(AminoAcidBase::getAminoAcidPtrByName("None") != nullptr);
        REQUIRE(AminoAcidBase::getAminoAcidPtrByName("Lysine") != nullptr);

        //expected to return null pointer
        REQUIRE(AminoAcidBase::getAminoAcidPtrByName(" ") == nullptr);
        REQUIRE(AminoAcidBase::getAminoAcidPtrByName("PREC") == nullptr);//this is ion name
    }

    SECTION("get amino acid pointer by letters"){
        //expected to return valid pointer
        AminoAcidPtr seleno = AminoAcidBase::getAminoAcidPtrByOneLetter("U");
        AminoAcidPtr none = AminoAcidBase::getAminoAcidPtrByOneLetter("X");
        
        AminoAcidPtr glutamic = AminoAcidBase::getAminoAcidPtrByThreeLetter("Glu");
        AminoAcidPtr phenyl = AminoAcidBase::getAminoAcidPtrByThreeLetter("Phe");

        //examine other information such as composition at that pointer address
        SECTION("examine amino acid pointer address"){
            REQUIRE(none->getAvgMass() == 0.0);
            REQUIRE(seleno->getComposition() == "C3H5NOSe");
            REQUIRE(glutamic->getMonoMass() == 129.042593088041);
            REQUIRE(seleno->getName() == "Selenocysteine");
            REQUIRE(phenyl->getOneLetter() == "F");
            REQUIRE(none->getThreeLetter() == "Xxx");
            REQUIRE(glutamic->getXmlElementName() == "amino_acid");
        }

        //expected to return null pointer
        REQUIRE(AminoAcidBase::getAminoAcidPtrByOneLetter("B") == nullptr);
        REQUIRE(AminoAcidBase::getAminoAcidPtrByOneLetter("4") == nullptr);

        REQUIRE(AminoAcidBase::getAminoAcidPtrByThreeLetter("Lol") == nullptr);
        REQUIRE(AminoAcidBase::getAminoAcidPtrByThreeLetter(" ") == nullptr);

    }

    SECTION("check if the list contains this amono acid"){
        REQUIRE(AminoAcidBase::containsName("Alanine"));
        REQUIRE(AminoAcidBase::containsOneLetter("K"));
        REQUIRE(AminoAcidBase::containsThreeLetter("Met"));

        REQUIRE(AminoAcidBase::containsName("Dopamine") == false);
        REQUIRE(AminoAcidBase::containsOneLetter(" ") == false);
        REQUIRE(AminoAcidBase::containsThreeLetter("Avd") == false);
    }

    SECTION ("get pointer from XML"){
        
    }
}

}


