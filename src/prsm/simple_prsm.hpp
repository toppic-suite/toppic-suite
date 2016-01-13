#ifndef PROT_PRSM_SIMPLE_PRSM_HPP_
#define PROT_PRSM_SIMPLE_PRSM_HPP_

#include <string>
#include <xercesc/dom/DOM.hpp>

#include "base/proteoform.hpp"
#include "base/xml_dom_document.hpp"
#include "spec/ms_header.hpp"

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
  double getScore(){return score_;}
  int getSpectrumId(){return spectrum_id_;}
  const std::string& getSpectrumScan(){return spectrum_scan_;}
  int getSpectrumNum() {return spectrum_num_;}
  double getPrecMass(){return prec_mass_;}
  void setPrecursorId(int precursorId) {precursor_id_ = precursorId;}
  int getPrecursorId(){return precursor_id_;}

  std::vector<double>& getNTruncShifts() {return n_trunc_shifts_;}
  std::vector<double>& getCTruncShifts() {return c_trunc_shifts_;}

  void setNTruncShifts(const std::vector<double> &shifts) {n_trunc_shifts_ = shifts;}
  void setCTruncShifts(const std::vector<double> &c_term_shifts);

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

  std::string seq_name_;
  std::string seq_desc_;
  double prot_mass_;

  double score_;

  std::vector<double> n_trunc_shifts_;
  std::vector<double> c_trunc_shifts_;
};

} /* namespace prot */

#endif /* SIMPLE_PRSM_HPP_ */
