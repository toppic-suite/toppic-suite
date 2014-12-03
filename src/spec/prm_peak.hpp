#ifndef PROT_PRM_PEAK_HPP_
#define PROT_PRM_PEAK_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_peak.hpp"
#include "spec/support_peak.hpp"
#include "spec/ms.hpp"
#include "spec/sp_para.hpp"
#include "deconv_ms.hpp"

namespace prot {

#define PRM_PEAK_TYPE_ORIGINAL 0
#define PRM_PEAK_TYPE_REVERSED 1

class PrmPeak : public Peak {
 public:
  PrmPeak(DeconvPeakPtr base_peak_ptr, int base_type, double mono_mass, 
          double score);

  void addNghbEdge(DeconvPeakPtr deconv_peak_ptr, double offset,
                   SPTypePtr peak_type, double score);

  int getNeighborSize(){return neighbor_list_.size();}

  DeconvPeakPtr getBasePeakPtr(){return base_peak_ptr_;}

  double getMonoMass(){return mono_mass_;}

  double getScore(){return score_;}

  double getStrictTolerance(){return strict_tolerance_;}

  int getBaseType(){return base_type_;}

  double getNStrictCRelaxTolerance(){return n_strict_c_relax_tolerance_;}

  double getNRelaxCStrictTolerance(){return n_relax_c_strict_tolerance_;}

  int getBreakType();

  void setStrictTolerance(double tolerance){strict_tolerance_ = tolerance;}

  void setNStrictCRelacTolerance(double tolerance){
    n_strict_c_relax_tolerance_ = tolerance;}

  void setNRelaxCStrictTolerance(double tolerance){
    n_relax_c_strict_tolerance_ = tolerance;}

 private:
  DeconvPeakPtr base_peak_ptr_;
  double mono_mass_;
  double score_;
  int base_type_;
  double strict_tolerance_;
  double n_strict_c_relax_tolerance_;
  double n_relax_c_strict_tolerance_;
  SupportPeakPtrVec neighbor_list_;
};

typedef std::shared_ptr<PrmPeak> PrmPeakPtr;
typedef std::vector<PrmPeakPtr> PrmPeakPtrVec;
typedef std::shared_ptr<Ms<PrmPeakPtr>> PrmMsPtr;

inline bool prmPeakUp(const PrmPeakPtr &a, const PrmPeakPtr &b){
  return a->getPosition() < b->getPosition();
}

PrmMsPtr createMsTwoPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, double new_prec_mass);

PrmMsPtr createMsSixPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, double new_prec_mass); 

PrmMsPtr createShiftMsSixPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, double new_prec_mass, 
                             double shift);

std::pair<std::vector<int>,std::vector<int>> getIntMassErrorList(PrmMsPtr prm_ms_ptr, double scale,
                                                                 bool n_strict, bool c_strict);

std::vector<double> getMassList(PrmMsPtr prm_ms_ptr);
std::vector<double> getScoreList(PrmMsPtr prm_ms_ptr);

} /* namespace prot */

#endif /* PRM_PEAK_HPP_ */
