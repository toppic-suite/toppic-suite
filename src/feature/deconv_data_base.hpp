#ifndef PROT_FEATURE_DECONV_DATA_BASE_HPP_
#define PROT_FEATURE_DECONV_DATA_BASE_HPP_

#include <memory>
#include <vector>

#include "spec/peak.hpp"
#include "feature/feature_mng.hpp"
#include "feature/deconv_data.hpp"

namespace prot {

class DeconvDataBase {
 public:
  static DeconvDataPtr getDataPtr(PeakPtrVec &peak_list, FeatureMngPtr
                                  mng_ptr); 
  static DeconvDataPtr getDataPtr(PeakPtrVec &peak_list, double max_mass, int
                                  max_charge, FeatureMngPtr mng_ptr);
};

}

#endif
