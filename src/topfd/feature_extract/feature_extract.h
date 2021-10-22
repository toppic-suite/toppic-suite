//
// Created by abbash on 10/8/2021.
//

#ifndef TOPPIC_FEATURE_EXTRACT_H
#define TOPPIC_FEATURE_EXTRACT_H

#include "ms/feature/frac_feature.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/sample_feature.hpp"
#include "ms/feature/sample_feature_writer.hpp"
#include "ms/spec/deconv_peak.hpp"

#include "topfd/feature_extract/python_caller.h"
#include "topfd/feature_extract/spectral_data_processor.h"
#include "topfd/feature_extract/mzml_data_processor.h"
#include "topfd/feature_extract/xic_generator.h"

namespace toppic {
namespace feature_extract {
int process(TopfdParaPtr para_ptr, std::vector<std::string> spec_file_lst);
}
}

#endif //TOPPIC_FEATURE_EXTRACT_H
