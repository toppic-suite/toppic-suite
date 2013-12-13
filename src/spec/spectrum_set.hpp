/*
 * spectrum_set.hpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#ifndef PROT_SPECTRUM_SET_HPP_
#define PROT_SPECTRUM_SET_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_ms.hpp"
#include "spec/sp_para.hpp"
#include "spec/extend_peak.hpp"
#include "spec/prm_peak.hpp"

namespace prot {
class SpectrumSet {
public:
	SpectrumSet(DeconvMsPtr sp,double delta,SpParaPtr sp_para,double shift);
	DeconvMsPtr getDeconvMs(){return deconv_sp_;}
	ExtendPeakMs getSpThree(){return extend_ms_three_;}
	PrmPeakMS getSpTwo(){return prm_ms_two_;}
	PrmPeakMS getSpSix(){return prm_ms_six_;}
	PrmPeakMS getSpShiftSix(){return prm_ms_shift_six_;}
	double getDelta(){return delta_;}
private:
	DeconvMsPtr deconv_sp_;
	double delta_;
	PrmPeakMS prm_ms_two_;
	ExtendPeakMs extend_ms_three_;
	PrmPeakMS prm_ms_six_;
	PrmPeakMS prm_ms_shift_six_;
};

typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

SpectrumSetPtr getSpectrumSet(DeconvMsPtr spectrum,double delta,SpParaPtr sp_para,double shift);

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
