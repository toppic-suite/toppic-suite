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
#include "base/base_data.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/extend_sp_para.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class PtmFastFilterMng {
 public:
  PtmFastFilterMng(std::string config_file_name) {
    base_data = BaseDataPtr(new BaseData(config_file_name));
  }

  PtmFastFilterMng(std::string config_file_name,
                   std::map<std::string, std::string> conf) {
    base_data = BaseDataPtr(new BaseData(config_file_name));
    spectrum_file_name_ = conf["spectrumFileName"];
    search_db_file_name_ = conf["databaseFileName"];
    activation_ptr_ = ActivationFactory::getBaseActivationPtrByName(
        conf["activation"]);
  }
  //Candidate protein number for each spectrum
  int ptm_fast_filter_result_num_ = 20;
  int db_block_size_ = 5000000;
  int ptm_fast_filter_scale_ = 100;
  //spectrum parameters
  double ppo_ = 0.000015;
  bool use_min_tolerance_ = true;
  double min_tolerance_ = 0.01;
  PeakTolerancePtr peak_tolerance_ = PeakTolerancePtr(
      new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_));

  int min_peak_num = 10;
  double min_mass = 50.0;

  // extend sp parameter
  double IM_ = MassConstant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list
  std::vector<double> ext_offsets_ { { 0, -IM_, IM_ } };
  double extend_min_mass_ = 5000;
  ExtendSpParaPtr extend_sp_para_ = ExtendSpParaPtr(
      new ExtendSpPara(extend_min_mass_, ext_offsets_));
  ActivationPtr activation_ptr_;

  SpParaPtr sp_para_ = SpParaPtr(
      new SpPara(min_peak_num, min_mass, peak_tolerance_, extend_sp_para_,
                 activation_ptr_));

  std::string search_db_file_name_;
  std::string res_file_name;
  std::string spectrum_file_name_;
  std::string output_file_ext_;

  BaseDataPtr base_data;
};

typedef std::shared_ptr<PtmFastFilterMng> PtmFastFilterMngPtr;

} /* namespace tools */

#endif /* PTM_FASTFILTER_MNG_HPP_ */
