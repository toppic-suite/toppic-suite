#include <boost/algorithm/string.hpp>
#include "prsm/prsm_para.hpp"

namespace prot {

PrsmPara::PrsmPara(std::map<std::string, std::string> &arguments) {
  search_db_file_name_ = arguments["databaseFileName"];
  spec_file_name_ = arguments["spectrumFileName"];
  log_file_name_ = arguments["logFileName"];

  fix_mod_residue_list_ = FixResidueFactory::getFixResiduePtrVec(arguments["cysteineProtection"]);

  std::string prot_mod_str = arguments["allowProtMod"];
  std::vector<std::string> strs;
  boost::split(strs, prot_mod_str, boost::is_any_of(","));
  for (size_t i = 0; i < strs.size(); i++) {
    ProtModPtr ptr = ProtModFactory::getBaseProtModPtrByName(strs[i]);
    allow_prot_mod_list_.push_back(ptr);
  }

  std::string activation_name = arguments["activation"];
  ActivationPtr activation_ptr 
      = ActivationFactory::getBaseActivationPtrByName(activation_name);

  double ppo = atoi(arguments["errorTolerance"].c_str())*0.000001;
  bool use_min_tolerance = true;
  double min_tolerance = 0.01;
  PeakTolerancePtr peak_tolerance_ptr = PeakTolerancePtr(
      new PeakTolerance(ppo, use_min_tolerance, min_tolerance));

  /** extend sp parameter */
  double IM = MassConstant::getIsotopeMass();
  /** the set of offsets used to expand the monoisotopic mass list */
  std::vector<double> ext_offsets {{0, -IM, IM}};
  double extend_min_mass = 5000;

  int min_peak_num = 10;
  double min_mass = 50.0;

  sp_para_ptr_ = SpParaPtr(new SpPara(min_peak_num, min_mass, extend_min_mass,
                                      ext_offsets, peak_tolerance_ptr, activation_ptr));

}

} /* namespace prot */
