//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef PROT_PRSM_SIMPLE_PRSM_HPP_
#define PROT_PRSM_SIMPLE_PRSM_HPP_

#include <string>
#include <vector>

#include <xercesc/dom/DOM.hpp>

#include "seq/proteoform.hpp"
#include "xml/xml_dom_document.hpp"
#include "spec/ms_header.hpp"

namespace toppic {

class SimplePrsm;
typedef std::shared_ptr<SimplePrsm>   SimplePrsmPtr;
typedef std::vector<SimplePrsmPtr>    SimplePrsmPtrVec;
typedef std::vector<SimplePrsmPtrVec> SimplePrsmPtrVec2D;

class SimplePrsm {
 public:
  SimplePrsm(MsHeaderPtr header_ptr, int spectrum_num,
             ProteoformPtr proteo_ptr, int score);

  SimplePrsm(MsHeaderPtr header_ptr, int spectrum_num,
             const std::string & seq_name,
             const std::string & seq_desc,
             int score);

  explicit SimplePrsm(xercesc::DOMElement* element);

  std::string getFileName () {return file_name_;}

  void setFileName(const std::string & fname) {file_name_ = fname;}

  std::string getSeqName() {return seq_name_;}

  std::string getSeqDesc() {return seq_desc_;}

  double getScore() {return score_;}

  int getSpectrumId() {return spectrum_id_;}

  const std::string& getSpectrumScan() {return spectrum_scan_;}

  int getSpectrumNum() {return spectrum_num_;}

  double getPrecMass() {return prec_mass_;}

  void setPrecursorId(int precursorId) {precursor_id_ = precursorId;}

  int getPrecursorId() {return precursor_id_;}

  double getProteoformMass() {return prot_mass_;}

  std::vector<double>& getNTruncShifts() {return n_trunc_shifts_;}

  std::vector<double>& getCTruncShifts() {return c_trunc_shifts_;}

  void setNTruncShifts(const std::vector<double> &shifts) {n_trunc_shifts_ = shifts;}

  void setCTruncShifts(const std::vector<double> &c_term_shifts);

  xercesc::DOMElement* toXml(XmlDOMDocument* xml_doc);

  std::vector<std::string> toStrVec();

  static std::string getXmlElementName() {return "simple_prsm";}

  static bool cmpScoreDec(const SimplePrsmPtr a, const SimplePrsmPtr b) {
    if (a->getScore() == b->getScore()) {
      return a->getSeqName() < b->getSeqName();
    } else {
      return a->getScore() > b->getScore();
    }
  }

  static bool cmpIdInc(const SimplePrsmPtr a, const SimplePrsmPtr b) {
    if (a->getSpectrumId() == b->getSpectrumId()) {
      return a->getSeqName() < b->getSeqName();
    } else {
      return a->getSpectrumId() < b->getSpectrumId();
    }
  }

  static bool cmpIdIncScoreDec(const SimplePrsmPtr a, const SimplePrsmPtr b) {
    if (a->getSpectrumId() < b->getSpectrumId()) {
      return true;
    } else if (a->getSpectrumId() > b->getSpectrumId()) {
      return false;
    } else {
      if (a->getScore() == b->getScore()) {
        return a->getSeqName() < b->getSeqName();
      } else {
        return a->getScore() > b->getScore();
      }
    }
  }

  static bool cmpNameIncScoreDec(const SimplePrsmPtr a, const SimplePrsmPtr b) {
    if (a->getSeqName() < b->getSeqName()) {
      return true;
    } else if (a->getSeqName() > b->getSeqName()) {
      return false;
    } else {
      return a->getScore() > b->getScore();
    }
  }

 private:
  std::string file_name_;
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

} /* namespace toppic */

#endif /* SIMPLE_PRSM_HPP_ */
