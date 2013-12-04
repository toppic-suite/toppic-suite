/*
 * extend_peak.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#ifndef PROT_EXTEND_PEAK_HPP_
#define PROT_EXTEND_PEAK_HPP_

#include <vector>
#include <memory>
#include <algorithm>

#include "spec/deconv_peak.hpp"
#include "spec/ms.hpp"
#include "spec/sp_para.hpp"
#include "spec/deconv_ms.hpp"

namespace prot {

class ExtendPeak : public Peak{
public:
	ExtendPeak();
	ExtendPeak(DeconvPeakPtr base_peak,double mono_mass,double score);
	DeconvPeakPtr getBasePeak(){return base_peak_;}
	double getMonoMass(){return mono_mass_;}
	double getScore(){return score_;}
	double getOrigTolerance(){return orig_tolerance_;}
	double getReverseTolerance(){return reverse_tolerance_;}
	void setOrigTolerance(double orig_tolerance){orig_tolerance_ = orig_tolerance;}
	void setReverseTolerance(double reverse_tolerance){reverse_tolerance_ = reverse_tolerance;}
private:
	DeconvPeakPtr base_peak_;
	double mono_mass_;
	double score_;
	double orig_tolerance_;
	double reverse_tolerance_;
};

typedef std::shared_ptr<ExtendPeak> ExtendPeakPtr;
typedef std::vector<ExtendPeakPtr> ExtendPeakPtrVec;

inline bool extendpeak_up(const ExtendPeakPtr p,ExtendPeakPtr n){
  return p->getPosition() < n->getPosition();
}

Ms<ExtendPeakPtr> getMsThree(DeconvMsPtr deconv_ms,double delta,SpParaPtr sp_para);

} /* namespace prot */

#endif /* EXTEND_PEAK_HPP_ */
