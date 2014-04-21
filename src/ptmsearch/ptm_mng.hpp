/*
 * ptm_mng.hpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PTM_MNG_HPP_
#define PROT_PTM_MNG_HPP_

#include <map>
#include "spec/peak_tolerance.hpp"
#include "base/mass_constant.hpp"
#include "base/trunc.hpp"
#include "base/base_data.hpp"
#include "spec/extend_sp_para.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class PtmMng {
 public :
  PtmMng(const std::string &conf_file_name, 
         const std::string &search_db_file_name, 
         const std::string &spectrum_file_name,
         const std::string &input_file_ext,
         const std::string &output_file_ext) {
    search_db_file_name_ = search_db_file_name;
    spectrum_file_name_ = spectrum_file_name;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;

    base_data_  = BaseDataPtr (new BaseData(conf_file_name));
    peak_tolerance_ = PeakTolerancePtr(
        new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));
    extend_sp_para_ = ExtendSpParaPtr(new ExtendSpPara(extend_min_mass_, ext_offsets_));
    sp_para_ = SpParaPtr(new SpPara(min_peak_num_, min_mass_, peak_tolerance_, 
                                    extend_sp_para_, base_data_->getActivationPtr())); 
  }

  PtmMng(std::map<std::string, std::string> arguments){
    spectrum_file_name_ = arguments["spectrumFileName"];
    search_db_file_name_ = arguments["databaseFileName"];
    n_unknown_shift_=atoi(arguments["shiftNumber"].c_str());
    ppo_=atoi(arguments["errorTolerance"].c_str())*0.000001;

    base_data_ = BaseDataPtr(new BaseData(arguments));
    peak_tolerance_ = PeakTolerancePtr(
        new PeakTolerance(ppo_,use_min_tolerance_,min_tolerance_));

    extend_sp_para_ = ExtendSpParaPtr(
        new ExtendSpPara(extend_min_mass_,ext_offsets_));
    sp_para_ = SpParaPtr( new SpPara(min_peak_num_, min_mass_, 
                   peak_tolerance_, extend_sp_para_, activation_));
  }

  BaseDataPtr base_data_ ;

  double ppo_ = 0.000015;
  bool use_min_tolerance_ = true;
  double min_tolerance_ = 0.01;
  PeakTolerancePtr peak_tolerance_;

  int min_peak_num_ = 10;
  double min_mass_ = 50.0;
  double IM_ = MassConstant::getIsotopeMass();
  std::vector<double> ext_offsets_ = {0.0,-IM_,IM_};
  double extend_min_mass_ = 5000;
  ExtendSpParaPtr extend_sp_para_;
  ActivationPtr activation_ = nullptr;
  SpParaPtr sp_para_;

  int n_report_ = 1;
  int n_unknown_shift_ =2;
  int n_known_shift_ = 0;

  int ptm_fast_filter_scale_ = 100;
  int n_top_diagonals_ = 20;
  double min_double_gap=0.25;
  int min_diagonal_gap_ = (int)(ptm_fast_filter_scale_ * min_double_gap);
  double extend_diagonal_error_tolerance_ = 0.5;
  double test_term_match_error_tolerance_ = 0.05;
  double align_prefix_suffix_shift_thresh_ = 300;

  double align_min_gap = 0.5;
  double large_shift_thresh = 300;
  double large_shift_panelty = 0;

  double adjust_prec_step_width_ = 0.005;

  std::string search_db_file_name_;
  std::string spectrum_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<PtmMng> PtmMngPtr;

} /* namespace prot */

#endif /* PTM_MNG_HPP_ */
