/*
 * ptm_mng.hpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PTM_MNG_HPP_
#define PROT_PTM_MNG_HPP_

#include "spec/peak_tolerance.hpp"
#include "base/mass_constant.hpp"
#include "base/trunc.hpp"
#include "base/base_data.hpp"
#include "spec/extend_sp_para.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class PtmMng {
public :
	PtmMng(std::string config_file_name);
	BaseDataPtr base_data_ ;

	double ppo_ = 0.000015;
	bool use_min_tolerance_ = true;
	double min_tolerance_ = 0.01;
	PeakTolerancePtr peank_tolerance_ = PeakTolerancePtr(new PeakTolerance(ppo_,use_min_tolerance_,min_tolerance_));

	int min_peak_num_ = 10;
	double min_mass_ = 50.0;
	double IM_ = MassConstant::getIsotopeMass();
	std::vector<double> extoffsets = {0.0,-IM_,IM_};
	double extend_thresh_ = 5000;
	ExtendSpParaPtr extend_sp_para_ = ExtendSpParaPtr(new ExtendSpPara(extend_thresh_,extoffsets));
	ActivationPtr activation_ = nullptr;
	SpParaPtr sp_para_ = SpParaPtr(new SpPara(min_peak_num_,min_mass_,peank_tolerance_,extend_sp_para_,activation_));

	//todo::need be inited;
	ProtModPtrVec allow_prot_N_mods_;
	ProtModPtrVec allow_prot_C_mods_;
	TruncPtrVec allow_prot_N_truncs_;
	TruncPtrVec allow_prot_C_truncs_;
	PtmPtrVec allow_pep_N_mods_;
	PtmPtrVec allow_pep_C_mods_;

	int n_report_ = 1;
	int n_unknown_shift_ =2;
	int n_known_shift_ = 0;

	int ptm_fast_filter_scale_ = 100;
	int n_top_diagonals_ = 20;
	double min_double_gap=0.25;
	int min_diagonal_gap_ = (int)(ptm_fast_filter_scale_ * min_double_gap);
	double extend_diagonal_error_tolerance_ = 0.5;
	double test_term_mod_error_toerance_ = 0.0001;
	double prefix_suffix_shift_thesh_ = 300;

	double align_min_gap = 0.5;
	double large_shift_thresh = 300;
	double large_shift_panelty = 0;

	double adjust_prec_step_width_ = 0.005;

	std::string search_db_file_name_;
	std::string res_file_name_;
	std::string spectrum_file_name_;
	std::string input_file_ext_;
	std::string output_file_ext_;
};

typedef std::shared_ptr<PtmMng> PtmMngPtr;

} /* namespace prot */

#endif /* PTM_MNG_HPP_ */
