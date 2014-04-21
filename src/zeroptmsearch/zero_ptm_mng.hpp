#ifndef PROT_ZERO_PTM_MNG_HPP_
#define PROT_ZERO_PTM_MNG_HPP_

#include <string>
#include <array>
#include <map>
#include <algorithm>

#include "base/mass_constant.hpp"
#include "base/residue.hpp"
#include "base/base_data.hpp"
#include "base/prot_mod.hpp"
#include "base/activation.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/extend_sp_para.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class ZeroPtmMng {
 public:
  ZeroPtmMng(const std::string &conf_file_name, 
             const std::string &search_db_file_name, 
             const std::string &spectrum_file_name,
             const std::string &output_file_ext) {
    search_db_file_name_ = search_db_file_name;
    spectrum_file_name_ = spectrum_file_name;
    output_file_ext_ = output_file_ext;

    base_data_ptr_  = BaseDataPtr (new BaseData(conf_file_name));
    peak_tolerance_ptr_ = PeakTolerancePtr(
        new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));
    extend_sp_para_ptr_ = ExtendSpParaPtr(new ExtendSpPara(extend_min_mass_, ext_offsets_));
    sp_para_ptr_ = SpParaPtr(new SpPara(min_peak_num_, min_mass_, peak_tolerance_ptr_, 
                              extend_sp_para_ptr_, base_data_ptr_->getActivationPtr())); 
  }

  ZeroPtmMng(std::map<std::string,std::string> arguments) {
    search_db_file_name_ = arguments["databaseFileName"];
    spectrum_file_name_ = arguments["spectrumFileName"];
    ppo_ = atoi(arguments["errorTolerance"].c_str())*0.000001;

    base_data_ptr_ = BaseDataPtr(new BaseData(arguments));
    peak_tolerance_ptr_ = PeakTolerancePtr(
        new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));

    extend_sp_para_ptr_ = ExtendSpParaPtr(new ExtendSpPara(extend_min_mass_, ext_offsets_));
    sp_para_ptr_ = SpParaPtr(new SpPara(min_peak_num_, min_mass_, peak_tolerance_ptr_,
                                        extend_sp_para_ptr_, base_data_ptr_->getActivationPtr()));
  }

  BaseDataPtr base_data_ptr_;

  /** zero ptm fast filtering */
  int zero_ptm_filter_result_num_ = 20;
  /** number of reported PrSMs for each spectrum */
  int report_num_ = 1;

  /** spectrum parameters */
  double ppo_ = 0.000015;
  bool use_min_tolerance_ = true;
  double min_tolerance_ = 0.01;
  PeakTolerancePtr peak_tolerance_ptr_;

  /** extend sp parameter */
  double IM_ = MassConstant::getIsotopeMass();
  /** the set of offsets used to expand the monoisotopic mass list */
  std::vector<double> ext_offsets_ {{0, -IM_, IM_}};
  double extend_min_mass_ = 5000;
  ExtendSpParaPtr extend_sp_para_ptr_;

  int min_peak_num_ = 10;
  double min_mass_ = 50.0;

  SpParaPtr sp_para_ptr_;

  /** recalibration is used in ZeroPtmSlowMatch */
  bool   do_recalibration_ = false;
  double recal_ppo_ = 0.000015; // 15 ppm
  bool   ms_one_ms_two_same_recal_ = true;

  std::string search_db_file_name_;
  std::string spectrum_file_name_;

  std::string output_file_ext_;
};

typedef std::shared_ptr<ZeroPtmMng> ZeroPtmMngPtr;

} /* namespace_prot */

#endif
