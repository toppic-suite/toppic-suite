#include "feature/deconv_data_base.hpp"
#include "feature/deconv_util.hpp"
#include "feature/deconv_one_sp.hpp"

namespace prot {

DeconvOneSp::DeconvOneSp(FeatureMngPtr mng_ptr): 
    mng_ptr_(mng_ptr) {
	}

void DeconvOneSp::setData(PeakPtrVec &peak_list) {
  data_ptr_ = DeconvDataBase::getDataPtr(peak_list, mng_ptr_);
}

void DeconvOneSp::preprocess() {
  if (mng_ptr_->estimate_min_inte_) {
    PeakPtrVec peak_list = data_ptr_->getPeakList();
    std::vector<double> intes;
    for (size_t i = 0; i < peak_list.size(); i++) {
      intes.push_back(peak_list[i]->getIntensity());
      double min_inte = DeconvUtil::getBaseLine(intes);
      mng_ptr_->setMinInte(min_inte);
    }
  }
}

}
