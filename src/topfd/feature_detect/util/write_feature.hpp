//
// Created by abbash on 8/30/22.
//

#ifndef TOPPIC_WRITE_FEATURE_HPP
#define TOPPIC_WRITE_FEATURE_HPP

#include <string>
#include <fstream>
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/feature/feature.hpp"


namespace toppic {
namespace write_feature {
    void writeFeatures(const std::string &output_file_name, const std::vector<Feature> &features);
}
}


#endif //TOPPIC_WRITE_FEATURE_HPP
