/*
 * simple_prsm.hpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

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
    SimplePrsm(MsHeaderPtr header,ProteoformPtr seq,int score);
    SimplePrsm(xercesc::DOMElement* element);
    std::string getSeqName(){return seq_name_;}
    ProteoformPtr getSeq(){return seq_;}
    int getSeqId(){return seq_id_;}
    double getScore(){return score_;}
    int getSpectrumId(){return spectrum_id_;}
    std::string getSpectrumScan(){return spectrum_scan_;}
    double getPrecMass(){return prec_mass_;}
    void setPrecursorId(int precursorId){precursor_id_ = precursorId;}
    int getPrecursorId(){return precursor_id_;}
    int compareTo(SimplePrsmPtr simple_prsm_ptr);
    void findSeq(std::vector<ProteoformPtr> seqs);
    xercesc::DOMElement* toXml(XmlDOMDocument* xml_doc);
    bool isMatch(MsHeaderPtr header);

private:
    int spectrum_id_;
    std::string spectrum_scan_;
    int precursor_id_;
    double prec_mass_;
    ProteoformPtr seq_;
    int seq_id_;
    std::string seq_name_;
    double score_;
};

SimplePrsmPtrVec findSimplePrsms(SimplePrsmPtrVec simple_prsm,
                                 MsHeaderPtr header);

SimplePrsmPtrVec readSimplePrsm(std::string filename);

inline bool simplePrsmDown(const SimplePrsmPtr p,SimplePrsmPtr n){
  return p->getScore() > n->getScore();
}

} /* namespace prot */

#endif /* SIMPLE_PRSM_HPP_ */
