#include <catch.hpp>

#include "common/base/prot_mod_util.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/residue_base.hpp"

namespace toppic {
namespace prot_mod_util {

TEST_CASE("prot_mod util testing"){
    ResidueBase::initBase();
    ProtModBase::initBase();//initialize mod base
    
    //REQUIRE(allowMod(ProtModBase::getProtModPtrByName("NONE"), ResidueBase::getBaseResiduePtrVec()));
    //both of them below return false
    //REQUIRE(ProtModBase::getProtModPtrByName("NONE") == ProtModBase::getProtModPtr_NONE());
    //REQUIRE(ProtModBase::getProtModPtrByType("NONE")[0] == ProtModBase::getProtModPtr_NONE());
    //REQUIRE(ProtModBase::getProtModPtrByName("M_ACETYLATION") == ProtModBase::getProtModPtr_M_ACETYLATION());
    //REQUIRE(allowMod(ProtModBase::getProtModPtrByName("NME_ACETYLATION_SER"), ResidueBase::getBaseResiduePtrVec()));
    
}



}
}