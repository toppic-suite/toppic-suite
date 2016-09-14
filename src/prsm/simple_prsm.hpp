// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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

  static bool cmpNameIncScoreDec(const SimplePrsmPtr a,SimplePrsmPtr b) {
    if (a->getSeqName() < b->getSeqName()) {
      return true;
    }
    else if (a->getSeqName() > b->getSeqName()) {
      return false;
    }
    else {
      return a->getScore() > b->getScore();
    }
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
