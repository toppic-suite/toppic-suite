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