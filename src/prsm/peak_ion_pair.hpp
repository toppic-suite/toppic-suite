#ifndef PROT_PRSM_PEAK_ION_PAIR_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_HPP_

#include "base/xml_dom_document.hpp"
#include "spec/ms_header.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"

namespace prot {

class PeakIonPair;
typedef std::shared_ptr<PeakIonPair> PeakIonPairPtr;

class PeakIonPair {
 public:
  PeakIonPair(MsHeaderPtr ms_header_ptr, ExtendPeakPtr real_peak_ptr, 
              TheoPeakPtr theo_peak_ptr); 

  MsHeaderPtr getMsHeaderPtr() {return ms_header_ptr_;}

  ExtendPeakPtr getRealPeakPtr() {return real_peak_ptr_;}

  TheoPeakPtr getTheoPeakPtr() {return theo_peak_ptr_;}

  void setId(int id) {id_ = id;}

  void appendRealPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

  void appendTheoPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

  static bool cmpRealPeakPosInc(const PeakIonPairPtr &a, const PeakIonPairPtr &b);

 private:
  int id_;
  MsHeaderPtr ms_header_ptr_;
  ExtendPeakPtr real_peak_ptr_;
  TheoPeakPtr theo_peak_ptr_;
};
typedef std::vector<PeakIonPairPtr> PeakIonPairPtrVec;
typedef std::vector<PeakIonPairPtrVec> PeakIonPairPtrVec2D;

}

#endif
