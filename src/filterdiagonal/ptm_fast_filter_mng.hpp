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
#include <map>

#include "base/mass_constant.hpp"
#include "base/activation.hpp"
#include "base/residue.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class PtmFastFilterMng {
 public:

  PtmFastFilterMng(std::map<std::string, std::string> arguments) {
    spectrum_file_name_ = arguments["spectrumFileName"];
    search_db_file_name_ = arguments["databaseFileName"];

    fix_mod_residue_list_ = FixResidueFactory::getFixResiduePtrVec(arguments["cysteineProtection"]);

    activation_ptr_ = ActivationFactory::getBaseActivationPtrByName(
        arguments["activation"]);
    ppo_=atoi(arguments["errorTolerance"].c_str())*0.000001;

    peak_tolerance_ = PeakTolerancePtr(
        new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));

    sp_para_ = SpParaPtr(
        new SpPara(min_peak_num, min_mass, extend_min_mass_,
                   ext_offsets_, peak_tolerance_, activation_ptr_));
  }

  std::string search_db_file_name_;
  std::string spectrum_file_name_;

  ResiduePtrVec fix_mod_residue_list_;
  // if activation ptr is null, activation types in file are used 
  ActivationPtr activation_ptr_;

  //spectrum parameters
  double ppo_ = 0.000015;
  bool use_min_tolerance_ = true;
  double min_tolerance_ = 0.01;
  PeakTolerancePtr peak_tolerance_;

  int min_peak_num = 10;
  double min_mass = 50.0;

  // extend sp parameter
  double IM_ = MassConstant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list
  std::vector<double> ext_offsets_ { { 0, -IM_, IM_ } };
  double extend_min_mass_ = 5000;

  SpParaPtr sp_para_;

  /** parameters for fast filteration */

  //Candidate protein number for each spectrum
  int ptm_fast_filter_result_num_ = 20;
  int db_block_size_ = 5000000;
  int ptm_fast_filter_scale_ = 100;

  std::string output_file_ext_;

};

typedef std::shared_ptr<PtmFastFilterMng> PtmFastFilterMngPtr;

} /* namespace tools */

#endif /* PTM_FASTFILTER_MNG_HPP_ */
