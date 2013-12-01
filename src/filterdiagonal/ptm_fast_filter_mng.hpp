/*
 * ptm_fastfilter_mng.hpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PTM_FASTFILTER_MNG_HPP_
#define PROT_PTM_FASTFILTER_MNG_HPP_

#include <vector>
#include <memory>

#include "peak_tolerance.hpp"
#include "mass_constant.hpp"
#include "extend_sp_para.hpp"
#include "activation.hpp"

namespace prot {

class PtmFastFilterMng {
public:
	//Candidate protein number for each spectrum
	int ptm_fast_filter_result_num_ = 20;
	int db_block_size_ = 5000000;
	int ptm_fast_filter_scale_ = 100;
	//spectrum parameters
	double ppo_ = 0.000015;
	bool use_min_tolerance_ = true;
	double min_tolerance_ = 0.01;
	PeakTolerance peak_tolerance_ = PeakTolerance(ppo_,use_min_tolerance_,min_tolerance_);

	int min_peak_num = 10;
	double min_mass =50.0;

	// extend sp parameter
	double IM_ = MassConstant::getIsotopeMass();
	// the set of offsets used to expand the monoisotopic mass list
	std::vector<double> ext_offsets_ {{0, -IM_, IM_}};
	double extend_min_mass_ = 5000;
	ExtendSpPara extend_sp_para_ = ExtendSpPara(extend_min_mass_, ext_offsets_);
	ActivationPtr activation__ptr_;

	//todo:SpPara

	std::string search_db_file_name_;
	std::string res_file_name;
	std::string spectrum_file_name_;
	std::string output_file_ext_;
};

typedef std::shared_ptr<PtmFastFilterMng> PtmFastFilterMngPtr;

} /* namespace tools */

#endif /* PTM_FASTFILTER_MNG_HPP_ */
