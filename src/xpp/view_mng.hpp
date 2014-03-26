/*
 * view_mng.hpp
 *
 *  Created on: Mar 25, 2014
 *      Author: xunlikun
 */

#ifndef VIEW_MNG_HPP_
#define VIEW_MNG_HPP_

namespace prot {

class ViewMng {
 public:
  ViewMng(std::string conf_file_name):
        base_data_ptr_ (new BaseData(conf_file_name)),
        peak_tolerance_ptr_ (
            new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_)),
        extend_sp_para_ptr_ (new ExtendSpPara(extend_min_mass_, ext_offsets_)),
        sp_para_ptr_(new SpPara(min_peak_num_, min_mass_, peak_tolerance_ptr_,
                                extend_sp_para_ptr_,
                                base_data_ptr_->getActivationPtr()))
    {}

  BaseDataPtr base_data_ptr_;

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

};

typedef std::shared_ptr<ViewMng> ViewMngPtr;
} /* namespace prot */

#endif /* VIEW_MNG_HPP_ */
