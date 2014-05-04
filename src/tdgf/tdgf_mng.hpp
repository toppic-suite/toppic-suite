#ifndef PROT_TDGF_MNG_HPP_
#define PROT_TDGF_MNG_HPP_

#include <map>

#include "base/mass_constant.hpp"
#include "base/base_data.hpp"
#include "base/prot_mod.hpp"
#include "base/activation.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class TdgfMng {
 public:
  TdgfMng(std::map<std::string,std::string> arguments) {
    search_db_file_name_=arguments["databaseFileName"];
    spectrum_file_name_=arguments["spectrumFileName"];

    fix_mod_residue_list_ = FixResidueFactory::getFixResiduePtrVec(arguments["cysteineProtection"]);

    ProtModPtr ptr = ProtModFactory::getBaseProtModPtrByName("NONE");
    allow_prot_mod_list_.push_back(ptr);
    ptr = ProtModFactory::getBaseProtModPtrByName("NME");
    allow_prot_mod_list_.push_back(ptr);
    ptr = ProtModFactory::getBaseProtModPtrByName("NME_ACETYLATION");
    allow_prot_mod_list_.push_back(ptr);

    std::string activation_type = arguments["activation"];
    activation_ptr_ = ActivationFactory::getBaseActivationPtrByName(
        activation_type);

    ppo_=atoi(arguments["errorTolerance"].c_str())*0.000001;

    peak_tolerance_ptr_ = PeakTolerancePtr(
        new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));

    sp_para_ptr_ = SpParaPtr(new SpPara(min_peak_num_, min_mass_, extend_min_mass_,
                                        ext_offsets_, peak_tolerance_ptr_, activation_ptr_));
  }

  std::string search_db_file_name_;
  std::string spectrum_file_name_;

  ResiduePtrVec fix_mod_residue_list_;
  ProtModPtrVec allow_prot_mod_list_;
  // if activation ptr is null, activation types in file are used 
  ActivationPtr activation_ptr_;

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

  int min_peak_num_ = 10;
  double min_mass_ = 50.0;

  SpParaPtr sp_para_ptr_;

  /** parameters for tdmg */

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
