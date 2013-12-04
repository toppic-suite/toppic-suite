/*
 * prm_peak.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PRM_PEAK_HPP_
#define PROT_PRM_PEAK_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_peak.hpp"
#include "spec/support_peak.hpp"
#include "spec/support_peak_type.hpp"
#include "spec/prm_peak_type.hpp"

namespace prot {

class PrmPeak : public Peak {
public:
	PrmPeak(DeconvPeakPtr base_peak,PrmPeakTypePtr base_type,double mono_mass,double score);
	void addNghbEdge(DeconvPeakPtr peak,double offset,SupportPeakTypePtr peak_type,double score);
	int getNeighborSize(){return neighbor_list_.size();}
	DeconvPeakPtr getBasePeak(){return base_peak_;}
	double getMonoMass(){return mono_mass_;}
	double getScr(){return score_;}
	double getStrictTolerance(){return strict_tolerance_;}
	PrmPeakTypePtr getPrmPeakType(){return base_type_;}
	double getNStrictCrelacTolerance(){return n_strict_c_relax_tolerance_;}
	double getNRelaxCStrictTolerance(){return n_relax_c_strict_tolerance_;}
	int getBreakType(SupportPeakTypePtrVec support_peak_type_list);
	void setNStrictCrelacTolerance(double tolerance){n_strict_c_relax_tolerance_ = tolerance;}
	void setNRelaxCStrictTolerance(double tolerance){n_relax_c_strict_tolerance_ = tolerance;}

private:
	DeconvPeakPtr base_peak_;
	double mono_mass_;
	double score_;
	PrmPeakTypePtr base_type_;
	double strict_tolerance_;
	double n_strict_c_relax_tolerance_;
	double n_relax_c_strict_tolerance_;
	SupportPeakPtrVec neighbor_list_;
};

typedef std::shared_ptr<PrmPeak> PrmPeakPtr;
typedef std::vector<PrmPeakPtr> PrmPeakPtrVec;

} /* namespace prot */

#endif /* PRM_PEAK_HPP_ */
