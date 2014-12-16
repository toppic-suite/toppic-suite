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
  PrmPeak(int spec_id, DeconvPeakPtr base_peak_ptr, int base_type, 
          double mono_mass, double score);

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

  int getSpecId() {return spec_id_;}

  int getPeakId() {return peak_id_;}

  int getBreakType();

  void setStrictTolerance(double tolerance){strict_tolerance_ = tolerance;}

  void setNStrictCRelacTolerance(double tolerance){
    n_strict_c_relax_tolerance_ = tolerance;}

  void setNRelaxCStrictTolerance(double tolerance){
    n_relax_c_strict_tolerance_ = tolerance;}

  void setPeakId(int peak_id) {peak_id_ = peak_id;}

 private:
  DeconvPeakPtr base_peak_ptr_;
  int peak_id_;
  int spec_id_;
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
typedef std::vector<PrmPeakPtrVec> PrmPeakPtrVec2D;
typedef std::shared_ptr<Ms<PrmPeakPtr>> PrmMsPtr;
typedef std::vector<PrmMsPtr> PrmMsPtrVec;

inline bool prmPeakUp(const PrmPeakPtr &a, const PrmPeakPtr &b){
  return a->getPosition() < b->getPosition();
}

/*
PrmMsPtr createMsTwoPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, double new_prec_mass);

PrmMsPtr createMsSixPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, double new_prec_mass); 

PrmMsPtr createShiftMsSixPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, double new_prec_mass, 
                             double shift);
                             */

PrmMsPtrVec createMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, SpParaPtr sp_para_ptr,
                              double prec_mono_mass);

PrmMsPtrVec createMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, SpParaPtr sp_para_ptr,
                              double prec_mono_mass);

PrmMsPtrVec createShiftMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, SpParaPtr sp_para_ptr, 
                                   double prec_mono_mass, double shift);

std::vector<std::pair<int, int>> getIntMassErrorList(const PrmMsPtrVec &prm_ms_ptr_vec, 
                                                     PeakTolerancePtr tole_ptr,
                                                     double scale, bool n_strict, bool c_strict);

PrmPeakPtrVec getPrmPeakPtrs(const PrmMsPtrVec &prm_ms_ptr_vec, PeakTolerancePtr tole_ptr);

/*

std::vector<double> getMassList(PrmMsPtr prm_ms_ptr);
std::vector<double> getMassList(const PrmMsPtrVec &prm_ms_ptr_vec);

std::vector<double> getScoreList(PrmMsPtr prm_ms_ptr);
*/

} /* namespace prot */

#endif /* PRM_PEAK_HPP_ */
