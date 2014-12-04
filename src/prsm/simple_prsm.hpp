#ifndef PROT_SIMPLE_PRSM_HPP_
#define PROT_SIMPLE_PRSM_HPP_

#include <string>
#include <xercesc/dom/DOM.hpp>

#include "htslib/faidx.h"

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
    SimplePrsm(MsHeaderPtr header_ptr,int spectrum_num,
               ProteoformPtr proteo_ptr, int score);
    SimplePrsm(xercesc::DOMElement* element);
    std::string getSeqName(){return seq_name_;}
    ProteoformPtr getProteoformPtr(){return proteo_ptr_;}
    ProteoformPtrVec getModProteoformPtrs() {return mod_proteo_ptrs_;}
    int getSeqId(){return seq_id_;}
    double getScore(){return score_;}
    int getSpectrumId(){return spectrum_id_;}
    const std::string& getSpectrumScan(){return spectrum_scan_;}
    int getSpectrumNum() {return spectrum_num_;}
    double getPrecMass(){return prec_mass_;}
    void setPrecursorId(int precursorId) {precursor_id_ = precursorId;}
    int getPrecursorId(){return precursor_id_;}
    xercesc::DOMElement* toXml(XmlDOMDocument* xml_doc);

    //to study
    bool isSameSpectrum(MsHeaderPtr header_ptr);
    bool isLargerSpectrumId(MsHeaderPtr header_ptr);
    void assignProteoformPtr(const ProteoformPtrVec &proteo_ptrs,
                             const ProteoformPtrVec2D &mod_proteo_2d_ptrs);

    void addProteoformPtr(faidx_t *fai, const ResiduePtrVec &residue_list,
                          const ProtModPtrVec &prot_mods);


private:
    int spectrum_id_;
    std::string spectrum_scan_;
    int precursor_id_;
    int spectrum_num_;
    double prec_mass_;

    ProteoformPtr proteo_ptr_;
    ProteoformPtrVec mod_proteo_ptrs_;
    int seq_id_;
    std::string seq_name_;
    std::string seq_desc_;
    double score_;
};

SimplePrsmPtrVec getMatchedSimplePrsms(const SimplePrsmPtrVec &simple_prsm_ptrs,
                                       MsHeaderPtr header);

SimplePrsmPtrVec readSimplePrsms(const std::string &filename);

SimplePrsmPtrVec getUniqueMatches(SimplePrsmPtrVec &match_ptrs);


inline bool simplePrsmDown(const SimplePrsmPtr a,SimplePrsmPtr b) {
  return a->getScore() > b->getScore();
}

inline bool simplePrsmSeqIdUpScoreDown(const SimplePrsmPtr a,SimplePrsmPtr b) {
  if (a->getSeqId() != b->getSeqId()) {
    return a->getSeqId() < b->getSeqId();
  }
  else {
    return a->getScore() > b->getScore();
  }
}

} /* namespace prot */

#endif /* SIMPLE_PRSM_HPP_ */
