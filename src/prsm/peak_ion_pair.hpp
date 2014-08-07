#ifndef PROT_PEAK_ION_PAIR_HPP_
#define PROT_PEAK_ION_PAIR_HPP_

#include "base/xml_dom_document.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"

namespace prot {

#define N_TERM_COVERAGE 0
#define C_TERM_COVERAGE 1
#define BOTH_TERM_COVERAGE 2

class PeakIonPair {
 public:
  PeakIonPair(ExtendPeakPtr real_peak_ptr, TheoPeakPtr theo_peak_ptr); 

  ExtendPeakPtr getRealPeakPtr() {return real_peak_ptr_;}

  TheoPeakPtr getTheoPeakPtr() {return theo_peak_ptr_;}

  void setId(int id) {id_ = id;}

  void appendRealPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

  void appendTheoPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

 private:
  int id_;
  ExtendPeakPtr real_peak_ptr_;
  TheoPeakPtr theo_peak_ptr_;
};

typedef std::shared_ptr<PeakIonPair> PeakIonPairPtr;
typedef std::vector<PeakIonPairPtr> PeakIonPairPtrVec;
typedef std::vector<PeakIonPairPtrVec> PeakIonPairPtrVec2D;

inline bool peakIonPairUp(const PeakIonPairPtr &a, const PeakIonPairPtr &b) {
  return a->getRealPeakPtr()->getBasePeakPtr()->getPosition() 
      < b->getRealPeakPtr()->getBasePeakPtr()->getPosition();
}

PeakIonPairPtrVec getMatchedPairs(const PeakIonPairPtrVec &pair_ptrs, 
                                  int peak_id);

PeakIonPairPtrVec getPeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                  const ExtendMsPtr &ms_three_ptr, 
                                  double min_mass);

double computePairConverage(const PeakIonPairPtrVec &pair_ptrs, int begin, 
                            int end, int coverage_type);

}

#endif
