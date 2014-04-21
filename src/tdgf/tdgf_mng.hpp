#ifndef PROT_TDGF_MNG_HPP_
#define PROT_TDGF_MNG_HPP_

#include <map>

#include "base/mass_constant.hpp"
#include "base/base_data.hpp"
#include "base/prot_mod.hpp"
#include "base/activation.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/extend_sp_para.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class TdgfMng {
 public:
  TdgfMng(const std::string &conf_file_name, 
          const std::string &search_db_file_name, 
          const std::string &spectrum_file_name,
          const std::string &input_file_ext,
          const std::string &output_file_ext) {
    search_db_file_name_ = search_db_file_name;
    spectrum_file_name_ = spectrum_file_name;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;

    base_data_ptr_  = BaseDataPtr (new BaseData(conf_file_name));
    peak_tolerance_ptr_ = PeakTolerancePtr(
        new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));
    extend_sp_para_ptr_ = ExtendSpParaPtr(new ExtendSpPara(extend_min_mass_, ext_offsets_));
    sp_para_ptr_ = SpParaPtr(new SpPara(min_peak_num_, min_mass_, peak_tolerance_ptr_, 
                                        extend_sp_para_ptr_, base_data_ptr_->getActivationPtr())); 
  }

  TdgfMng(std::map<std::string,std::string> arguments) {
    search_db_file_name_=arguments["databaseFileName"];
    spectrum_file_name_=arguments["spectrumFileName"];
    ppo_=atoi(arguments["errorTolerance"].c_str())*0.000001;

    base_data_ptr_ = BaseDataPtr (new BaseData(arguments)),
    peak_tolerance_ptr_ = PeakTolerancePtr(
        new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));
    extend_sp_para_ptr_ = ExtendSpParaPtr(new ExtendSpPara(extend_min_mass_, ext_offsets_));
    sp_para_ptr_ = SpParaPtr(new SpPara(min_peak_num_, min_mass_, peak_tolerance_ptr_, 
                              extend_sp_para_ptr_, base_data_ptr_->getActivationPtr())); 
  }

  BaseDataPtr base_data_ptr_;

  std::string search_db_file_name_;
  std::string spectrum_file_name_;

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

  /** PrSM filter */
  int comp_evalue_min_peak_num_ = 4;

  /** dp table */
  // number of mass shift
  int unexpected_shift_num_ = 2;
  double double_to_int_constant_ = 274.335215;
  double max_sp_prec_mass_ = 51000;
  int max_table_height_ = 128;
  int min_height_ = 10;

  /* Semi alignment type determination */
  /* a prsm with a shift < 300 at n-terminus is treated as a prefix */
  double prefix_suffix_shift_thresh_ = 300;

  /**
   * Postprocessing: adjustment makes it more conservative to identify PrSMs
   * with multiple shifts
   */
  // the following three values should be adjusted
  double multiple_shift_adjustment_ = 4;
  double multiple_shfit_adjustment_base_ = 10;
  double min_adjustment_ = 100;

  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<TdgfMng> TdgfMngPtr;

}
#endif
