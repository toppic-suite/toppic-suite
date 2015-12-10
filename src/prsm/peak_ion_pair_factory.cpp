#include <cmath>

#include "base/logger.hpp"
#include "base/algorithm.hpp"
#include "spec/extend_ms.hpp"
#include "spec/theo_peak.hpp"
#include "spec/theo_peak_factory.hpp"
#include "spec/theo_peak_util.hpp"
#include "prsm/peak_ion_pair_factory.hpp"

namespace prot {

void findPairs(ExtendMsPtr ms_three_ptr, TheoPeakPtrVec &theo_peak_ptrs, 
               int bgn, int end, PeakIonPairPtrVec &pair_ptrs) {
  std::sort(theo_peak_ptrs.begin(), theo_peak_ptrs.end(), TheoPeak::cmpPosInc);
  std::vector<double> ms_masses = ExtendMs::getExtendMassVec(ms_three_ptr);
  std::vector<double> theo_masses = TheoPeakUtil::getTheoMassVec(theo_peak_ptrs);

  size_t i = 0;
  size_t j = 0;
  while (i < ms_masses.size() && j < theo_masses.size()) {
    double deviation = ms_masses[i] - theo_masses[j];
    IonPtr ion_ptr = theo_peak_ptrs[j]->getIonPtr();
    double err = ms_three_ptr->getPeakPtr(i)->getOrigTolerance();
    if (ion_ptr->getPos() >= bgn && ion_ptr->getPos() <= end) {
      if (std::abs(deviation) <= err) {
        PeakIonPairPtr pair_ptr = PeakIonPairPtr(new PeakIonPair(
                ms_three_ptr->getMsHeaderPtr(), ms_three_ptr->getPeakPtr(i), theo_peak_ptrs[j]));
        pair_ptrs.push_back(pair_ptr);
      }
    }
    if (increaseIJ(i, j, deviation, err, ms_masses, theo_masses)) {
      i++;
    } else {
      j++;
    }
  }
}

/* parameter min_mass is necessary */
PeakIonPairPtrVec PeakIonPairFactory::genePeakIonPairs (const ProteoformPtr &proteoform_ptr, 
                                                        const ExtendMsPtr &ms_three_ptr, 
                                                        double min_mass) {
  ActivationPtr activation_ptr 
      = ms_three_ptr->getMsHeaderPtr()->getActivationPtr();

  TheoPeakPtrVec theo_peaks = TheoPeakFactory::geneProteoformTheoPeak(proteoform_ptr, 
                                                                      activation_ptr, 
                                                                      min_mass);

  PeakIonPairPtrVec pair_ptrs;
  findPairs(ms_three_ptr, theo_peaks, 0, proteoform_ptr->getLen(), pair_ptrs);
  return pair_ptrs;

}

PeakIonPairPtrVec PeakIonPairFactory::genePeakIonPairs(const ProteoformPtr &proteoform_ptr,
                                                       const ExtendMsPtrVec &ms_ptr_vec, 
                                                       double min_mass) {

  PeakIonPairPtrVec pair_ptrs;
  for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
    PeakIonPairPtrVec pair_ptr_tmp = genePeakIonPairs(proteoform_ptr, ms_ptr_vec[i],
                                                      min_mass);
    pair_ptrs.insert(pair_ptrs.end(), pair_ptr_tmp.begin(), pair_ptr_tmp.end());
  }
  return pair_ptrs;
}

}
