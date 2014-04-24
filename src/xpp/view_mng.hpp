/*
 * view_mng.hpp
 *
 *  Created on: Mar 25, 2014
 *      Author: xunlikun
 */

#ifndef VIEW_MNG_HPP_
#define VIEW_MNG_HPP_

#include <string>
#include <algorithm>
#include <map>

namespace prot {

class ViewMng {
 public:
  ViewMng(std::map<std::string,std::string> arguments) {
    spectrum_file_ = arguments["spectrumFileName"];
    database_file_=arguments["databaseFileName"];

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


    ppo_ = atoi(arguments["errorTolerance"].c_str())*0.000001;
    peak_tolerance_ptr_ = PeakTolerancePtr(new PeakTolerance(ppo_,use_min_tolerance_, min_tolerance_));
    extend_sp_para_ptr_=ExtendSpParaPtr(new ExtendSpPara(extend_min_mass_, ext_offsets_));
    sp_para_ptr_= SpParaPtr(new SpPara(min_peak_num_, min_mass_, peak_tolerance_ptr_,
                                       extend_sp_para_ptr_, activation_ptr_));
  }
  //std::map<std::string,std::string> arguments_;
  std::string spectrum_file_;
  std::string database_file_;

  ResiduePtrVec fix_mod_residue_list_;
  ProtModPtrVec allow_prot_mod_list_;
  // if activation ptr is null, activation types in file are used 
  ActivationPtr activation_ptr_;

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

  std::string html_path_="html/";
  std::string xml_path_="xml/";

};

typedef std::shared_ptr<ViewMng> ViewMngPtr;
} /* namespace prot */

#endif /* VIEW_MNG_HPP_ */
