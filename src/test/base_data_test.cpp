#include <catch.hpp>

#include "common/base/base_data.hpp"

namespace toppic {

namespace base_data {

TEST_CASE("check if all base data are initialized", "[classic]"){
    init();
    REQUIRE(base_data_init_);
}

}

}

