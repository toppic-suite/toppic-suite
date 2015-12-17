#ifndef PROT_SPEC_PRM_PEAK_HPP_
#define PROT_SPEC_PRM_PEAK_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_peak.hpp"
#include "spec/support_peak.hpp"
#include "spec/prm_base_type.hpp"
#include "spec/prm_break_type.hpp"

namespace prot {

class PrmPeak;
typedef std::shared_ptr<PrmPeak> PrmPeakPtr;

class PrmPeak : public Peak {
 public:
  PrmPeak(int spec_id, DeconvPeakPtr base_peak_ptr,
          PrmBaseTypePtr base_type, 
          double mono_mass, double score);

  void addNghbEdge(DeconvPeakPtr deconv_peak_ptr, double offset,
                   SPTypePtr peak_type, double score);

  int getNeighborSize(){return neighbor_list_.size();}

  DeconvPeakPtr getBasePeakPtr(){return base_peak_ptr_;}

  double getMonoMass(){return mono_mass_;}

  double getScore(){return score_;}

  double getStrictTolerance(){return strict_tolerance_;}

  PrmBaseTypePtr getBaseTypePtr(){return base_type_;}

  double getNStrictCRelaxTolerance(){return n_strict_c_relax_tolerance_;}

  double getNRelaxCStrictTolerance(){return n_relax_c_strict_tolerance_;}

  int getSpecId() {return spec_id_;}

  int getPeakId() {return peak_id_;}

  PrmBreakTypePtr getBreakType();

  void setStrictTolerance(double tolerance){strict_tolerance_ = tolerance;}

  void setNStrictCRelacTolerance(double tolerance){
    n_strict_c_relax_tolerance_ = tolerance;}

  void setNRelaxCStrictTolerance(double tolerance){
    n_relax_c_strict_tolerance_ = tolerance;}

  void setPeakId(int peak_id) {peak_id_ = peak_id;}

  static bool cmpPosIncrease(const PrmPeakPtr &a, const PrmPeakPtr &b){
    return a->getPosition() < b->getPosition();
  }

 private:
  int spec_id_;
  DeconvPeakPtr base_peak_ptr_;
  PrmBaseTypePtr base_type_;
  int peak_id_;
  double mono_mass_;
  double score_;
  double strict_tolerance_;
  double n_strict_c_relax_tolerance_;
  double n_relax_c_strict_tolerance_;
  SupportPeakPtrVec neighbor_list_;
};

typedef std::vector<PrmPeakPtr> PrmPeakPtrVec;
typedef std::vector<PrmPeakPtrVec> PrmPeakPtrVec2D;


} /* namespace prot */

#endif /* PRM_PEAK_HPP_ */
