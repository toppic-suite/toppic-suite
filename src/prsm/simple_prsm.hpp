#ifndef PROT_PRSM_SIMPLE_PRSM_HPP_
#define PROT_PRSM_SIMPLE_PRSM_HPP_

#include <string>
#include <xercesc/dom/DOM.hpp>

#include "base/proteoform.hpp"
#include "base/xml_dom_document.hpp"
#include "spec/ms_header.hpp"
#include "prsm/n_term_shift.hpp"

namespace prot {

class SimplePrsm;
typedef std::shared_ptr<SimplePrsm> SimplePrsmPtr;
typedef std::vector<SimplePrsmPtr> SimplePrsmPtrVec;
typedef std::vector<SimplePrsmPtrVec> SimplePrsmPtrVec2D;

class SimplePrsm {
 public:
  SimplePrsm(MsHeaderPtr header_ptr,int spectrum_num,
             ProteoformPtr proteo_ptr, int score);
  SimplePrsm(xercesc::DOMElement* element);
  std::string getSeqName(){return seq_name_;}
  std::string getSeqDesc(){return seq_desc_;}
  //ProteoformPtr getProteoformPtr(){return proteo_ptr_;}
  //ProteoformPtrVec getModProteoformPtrs() {return mod_proteo_ptrs_;}
  double getScore(){return score_;}
  int getSpectrumId(){return spectrum_id_;}
  const std::string& getSpectrumScan(){return spectrum_scan_;}
  int getSpectrumNum() {return spectrum_num_;}
  double getPrecMass(){return prec_mass_;}
  void setPrecursorId(int precursorId) {precursor_id_ = precursorId;}
  int getPrecursorId(){return precursor_id_;}
  NTermShiftPtrVec getNTermShiftPtrVec() {return n_term_shifts_;}

  xercesc::DOMElement* toXml(XmlDOMDocument* xml_doc);

  static std::string getXmlElementName() {return "simple_prsm";}

  static bool cmpScoreDec(const SimplePrsmPtr a,SimplePrsmPtr b) {
    return a->getScore() > b->getScore();
  }

  static bool cmpNameInc(const SimplePrsmPtr a,SimplePrsmPtr b) {
    return a->getSeqName() < b->getSeqName();
  }

 private:
  int spectrum_id_;
  std::string spectrum_scan_;
  int precursor_id_;
  int spectrum_num_;
  double prec_mass_;

  //ProteoformPtr proteo_ptr_;
  //ProteoformPtrVec mod_proteo_ptrs_;
  std::string seq_name_;
  std::string seq_desc_;
  double score_;

  NTermShiftPtrVec n_term_shifts_;
};

} /* namespace prot */

#endif /* SIMPLE_PRSM_HPP_ */
