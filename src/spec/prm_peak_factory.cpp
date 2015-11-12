
#include "base/ion_type_base.hpp"
#include "spec/prm_peak_factory.hpp"

namespace prot {

PrmPeakPtr PrmPeakFactory::getZeroPeakPtr(int spec_id, double prec_mono_mass, 
                                          PeakTolerancePtr tole_ptr, double score) {
  //zero_peak
  DeconvPeakPtr zero_peak_ptr = DeconvPeakPtr(new DeconvPeak(-1,0,0,0));
  PrmPeakPtr prm_peak_ptr(new PrmPeak(spec_id, zero_peak_ptr,
                                      PrmBaseType::ORIGINAL,0, score));
  // set tolerance 
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  return prm_peak_ptr;
}

PrmPeakPtr PrmPeakFactory::getPrecPeakPtr(int spec_id, double prec_mono_mass, 
                                          PeakTolerancePtr tole_ptr, double score) {
  //prec_peak
  double prec_peak_shift = IonTypeBase::getIonTypePtr_PREC()->getShift();
  double prec_peak_mass = prec_mono_mass -
      prec_peak_shift;
  DeconvPeakPtr prec_peak_ptr(new DeconvPeak(-1,prec_peak_mass, 0,0));
  PrmPeakPtr prm_peak_ptr = PrmPeakPtr(
      new PrmPeak(spec_id,prec_peak_ptr,
                  PrmBaseType::ORIGINAL,prec_peak_mass, score));
  // set tolerance
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(0));
  return prm_peak_ptr;
}

} /* namespace prot */
