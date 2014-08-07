#ifndef PROT_SIMPLE_PRSM_HPP_
#define PROT_SIMPLE_PRSM_HPP_

#include <string>
#include <xercesc/dom/DOM.hpp>

#include "base/proteoform.hpp"
#include "spec/ms_header.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class SimplePrsm;
typedef std::shared_ptr<SimplePrsm> SimplePrsmPtr;
typedef std::vector<SimplePrsmPtr> SimplePrsmPtrVec;
typedef std::vector<SimplePrsmPtrVec> SimplePrsmPtrVec2D;

class SimplePrsm {
public:
    SimplePrsm(MsHeaderPtr header_ptr,ProteoformPtr proteo_ptr, int score);
    SimplePrsm(xercesc::DOMElement* element);
    std::string getSeqName(){return seq_name_;}
    ProteoformPtr getProteoformPtr(){return proteo_ptr_;}
    int getSeqId(){return seq_id_;}
    double getScore(){return score_;}
    int getSpectrumId(){return spectrum_id_;}
    const std::string& getSpectrumScan(){return spectrum_scan_;}
    double getPrecMass(){return prec_mass_;}
    void setPrecursorId(int precursorId) {precursor_id_ = precursorId;}
    int getPrecursorId(){return precursor_id_;}
    xercesc::DOMElement* toXml(XmlDOMDocument* xml_doc);

    //to study
    bool isSameSpectrum(MsHeaderPtr header_ptr);
    bool isLargerSpectrumId(MsHeaderPtr header_ptr);
    void assignProteoformPtr(const std::vector<ProteoformPtr> &proteo_ptrs);

private:
    int spectrum_id_;
    std::string spectrum_scan_;
    int precursor_id_;
    double prec_mass_;

    ProteoformPtr proteo_ptr_;
    int seq_id_;
    std::string seq_name_;
    double score_;
};

SimplePrsmPtrVec getMatchedSimplePrsms(const SimplePrsmPtrVec &simple_prsm_ptrs,
                                       MsHeaderPtr header);

SimplePrsmPtrVec readSimplePrsms(const std::string &filename);

inline bool simplePrsmDown(const SimplePrsmPtr a,SimplePrsmPtr b) {
  return a->getScore() > b->getScore();
}

} /* namespace prot */

#endif /* SIMPLE_PRSM_HPP_ */
