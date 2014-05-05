#include "base/logger.hpp"
#include "base/algorithm.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

void PeakIonPair::appendPeakToXml(XmlDOMDocument* xml_doc, 
                                  xercesc::DOMElement* parent) {

  xercesc::DOMElement* element = xml_doc->createElement("matched_peak");
  std::string str = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName();
  xml_doc->addElement(element, "ion_name", str.c_str());
  str = convertToString(real_peak_ptr_->getBasePeakPtr()->getId());
  xml_doc->addElement(element, "peak_id", str.c_str());
  str = convertToString(real_peak_ptr_->getBasePeakPtr()->getCharge());
  xml_doc->addElement(element, "peak_charge", str.c_str());
  parent->appendChild(element);
}

void PeakIonPair::appendIonToXml(XmlDOMDocument* xml_doc, 
                                 xercesc::DOMElement* parent) {
  int pos=4;
  xercesc::DOMElement* element = xml_doc->createElement("matched_ion");
  std::string str 
      = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName().substr(0,1);
  xml_doc->addElement(element, "type", str.c_str());
  str = convertToString(theo_peak_ptr_->getShift());
  xml_doc->addElement(element, "match_shift", str.c_str()); 
  str = convertToString(real_peak_ptr_->getMonoMass(),pos);
  xml_doc->addElement(element, "adjusted_mass", str.c_str()); 
  str = convertToString(theo_peak_ptr_->getIonPtr()->getDisplayPos());
  xml_doc->addElement(element, "ion_display_position", str.c_str());
  str = convertToString(theo_peak_ptr_->getIonPtr()->getPos());
  xml_doc->addElement(element, "ion_left_position", str.c_str());
  double error = real_peak_ptr_->getMonoMass() - theo_peak_ptr_->getModMass();
  str = convertToString(error,pos);
  xml_doc->addElement(element, "mass_error", str.c_str()); 
  str = convertToString(error * 1000000 / real_peak_ptr_->getMonoMass(),pos-2);
  xml_doc->addElement(element, "ppm", str.c_str()); 
  parent->appendChild(element);
}

PeakIonPairPtrVec getMatchedPairs(const PeakIonPairPtrVec &pairs, int peak_id) {
  PeakIonPairPtrVec selected_pairs;
  for (unsigned int i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getRealPeakPtr()->getBasePeakPtr()->getId() == peak_id) {
      selected_pairs.push_back(pairs[i]);
    }
  }
  return selected_pairs;
}

void findPairs(ExtendMsPtr ms_three_ptr, TheoPeakPtrVec &theo_peaks, 
               int bgn, int end, PeakIonPairPtrVec &pairs) {
  std::sort(theo_peaks.begin(), theo_peaks.end(), theoPeakUp);
  std::vector<double> ms_masses = getExtendMassVec(ms_three_ptr);
  std::vector<double> theo_masses = getTheoMassVec(theo_peaks);

  unsigned int i = 0;
  unsigned int j = 0;
  while (i < ms_masses.size() && j < theo_masses.size()) {
    double deviation = ms_masses[i] - theo_masses[j];
    IonPtr ion_ptr = theo_peaks[j]->getIonPtr();
    double err = ms_three_ptr->getPeakPtr(i)->getOrigTolerance();
    if (ion_ptr->getPos() >= bgn && ion_ptr->getPos() <= end) {
      if (std::abs(deviation) <= err) {
        PeakIonPairPtr pair_ptr 
            = PeakIonPairPtr(new PeakIonPair(ms_three_ptr->getPeakPtr(i), 
                                             theo_peaks[j]));
        pairs.push_back(pair_ptr);
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
PeakIonPairPtrVec getPeakIonPairs (const ProteoformPtr &proteoform_ptr, 
                                   const ExtendMsPtr &ms_three_ptr, 
                                   double min_mass) {
  ActivationPtr activation_ptr 
      = ms_three_ptr->getHeaderPtr()->getActivationPtr();

  TheoPeakPtrVec theo_peaks = getProteoformTheoPeak(proteoform_ptr, 
                                                    activation_ptr, 
                                                    min_mass);

  PeakIonPairPtrVec pairs;
  findPairs(ms_three_ptr, theo_peaks, 0, proteoform_ptr->getLen(), pairs);
  return pairs;

}

double computePairConverage(const PeakIonPairPtrVec &pairs, int begin, 
                            int end, int coverage_type) {
  int total_num = end - begin  + 1;
  if (total_num <= 0) {
    return 0.0;
  }
  std::vector<bool> is_cov(total_num);
  for (unsigned int i  = 0; i < pairs.size(); i++) {
    IonPtr ion_ptr = pairs[i]->getTheoPeakPtr()->getIonPtr();
    bool cov = false;
    if (coverage_type == N_TERM_COVERAGE) {
      if (ion_ptr->getIonTypePtr()->isNTerm()) {
        cov = true;
      }
    }
    else if (coverage_type == C_TERM_COVERAGE) {
      if (!ion_ptr->getIonTypePtr()->isNTerm()) {
        cov = true;
      }
    }
    else if (coverage_type = BOTH_TERM_COVERAGE) {
      cov = true;
    }
    if (cov) {
      int pos = ion_ptr->getPos();
      if (pos >= begin && pos <= end) {
        is_cov[pos - begin] = true;
      }
    }
  }
  int cov_num = 0;
  for (unsigned int i = 0; i < is_cov.size(); i++) {
    if (is_cov[i]) {
      cov_num++;
    }
  }
  return cov_num/(double)total_num;

}

}
