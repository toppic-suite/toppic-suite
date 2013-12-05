#ifndef PROT_ZERO_PTM_MNG_HPP_
#define PROT_ZERO_PTM_MNG_HPP_

#include <string>
#include <array>
#include <map>

#include "base/mass_constant.hpp"
#include "base/residue.hpp"
#include "base/base_data.hpp"
#include "base/prot_mod.hpp"
#include "base/activation.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/extend_sp_para.hpp"

namespace prot {

class ZeroPtmMng {
 public:
  ZeroPtmMng(const char* conf_file_name);

  BaseDataPtr base_data_ptr_;

  /** zero ptm fast filtering */
  int zero_ptm_filter_result_num_ = 20;
  /** number of reported PrSMs for each spectrum */
  int report_num_ = 1;

  /** spectrum parameters */
  double ppo_ = 0.000015;
  bool use_min_tolerance_ = true;
  double min_tolerance_ = 0.01;
  PeakTolerance peak_tolerance_ = PeakTolerance (ppo_, use_min_tolerance_, min_tolerance_);

  /** extend sp parameter */
  double IM_ = MassConstant::getIsotopeMass();
  /** the set of offsets used to expand the monoisotopic mass list */
  std::vector<double> ext_offsets_ {{0, -IM_, IM_}};
  double extend_min_mass_ = 5000;
  ExtendSpPara extend_sp_para_ = ExtendSpPara(extend_min_mass_, ext_offsets_);

  int min_peak_num_ = 10;
  double min_mass_ = 50.0;

  ActivationPtr activation__ptr_;

  //SpPara spPara = new SpPara(minPeakNum, minMass, peakTolerance, extendSpPara, activationType);

  /** recalibration is used in ZeroPtmSlowMatch */
  bool   do_recalibration_ = true;
  double recal_ppo_ = 0.000015; // 15 ppm
  bool   ms_one_ms_two_same_recal_ = true;

  /** allowed N-terminal modifications */
  std::vector<ProtMod> allow_prot_mods_;
  std::map<std::string, ResiduePtr> fixed_ptms_;
  std::map<std::string, ResiduePtrVec> var_ptms_;

  std::string search_db_file_name_;
  std::string spectrum_file_name_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<ZeroPtmMng> ZeroPtmMngPtr;

void zeroPtmSearchProcess(ZeroPtmMngPtr);

} /* namespace_prot */

#endif
