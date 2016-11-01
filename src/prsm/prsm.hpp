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


#ifndef PROT_PRSM_PRSM_HPP_
#define PROT_PRSM_PRSM_HPP_

#include <string>

#include "base/extreme_value.hpp"
#include "base/proteoform.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/extend_ms.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class Prsm;
typedef std::shared_ptr<Prsm> PrsmPtr;

class Prsm {
 public:
  Prsm(ProteoformPtr proteoform_ptr, const DeconvMsPtrVec &deconv_ms_ptr_vec, 
       double adjusted_prec_mass, SpParaPtr sp_para_ptr);

  Prsm(xercesc::DOMElement* element, FastaIndexReaderPtr reader_ptr, 
       const ModPtrVec &fix_mod_list);

  int getPrsmId() {return prsm_id_;}

  int getSpectrumId() {return spectrum_id_;}

  std::string getSpectrumScan() {return spectrum_scan_;}

  int getPrecursorId() {return precursor_id_;}

  int getPrecFeatureId() {return prec_feature_id_;}
  
  double getPrecFeatureInte() {return prec_feature_inte_;}

  double getOriPrecMass() {return ori_prec_mass_;}

  double getAdjustedPrecMass() {return adjusted_prec_mass_;}

  ProteoformPtr getProteoformPtr() {return proteoform_ptr_;}

  ExtremeValuePtr getExtremeValuePtr() {return extreme_value_ptr_;} 

  double getFdr() {return fdr_;}

  double getProteoformFdr() {return proteoform_fdr_;}

  DeconvMsPtrVec getDeconvMsPtrVec() {return deconv_ms_ptr_vec_;}

  ExtendMsPtrVec getRefineMsPtrVec() {return refine_ms_three_vec_;}

  double getMatchPeakNum() {return match_peak_num_;}

  double getMatchFragNum() {return match_fragment_num_;}

  double getNormMatchFragNum();

  // ExtremeValue related functions
  double getEValue();

  double getPValue();

  double getOneProtProb();

  // set 
  void setPrsmId(int id) {prsm_id_ = id;}

  void setSpectrumId(int spectrum_id) {spectrum_id_ = spectrum_id;}

  void setSpectrumScan(std::string spectrum_scan) {spectrum_scan_ = spectrum_scan;}

  void setPrecurorId(int precursor_id) {precursor_id_ = precursor_id;}

  void setPrecFeatureId(int prec_feature_id) {prec_feature_id_ = prec_feature_id;}

  void setPrecFeatureInte(double prec_feature_inte) {prec_feature_inte_ = prec_feature_inte_;}

  void setOriPrecMass(double prec_mass) {ori_prec_mass_ = prec_mass;}

  void setProteoformPtr(ProteoformPtr proteoform) {proteoform_ptr_=proteoform;}

  void setExtremeValuePtr(ExtremeValuePtr ev_ptr) {extreme_value_ptr_ = ev_ptr;}

  void setFdr(double fdr) {fdr_ = fdr;}

  void setProteoformFdr(double proteoform_fdr) {proteoform_fdr_ = proteoform_fdr;}

  void setDeconvMsPtrVec(DeconvMsPtrVec ms_vec) {deconv_ms_ptr_vec_=ms_vec;}

  void setRefineMsVec(ExtendMsPtrVec refine_ms_three_vec){refine_ms_three_vec_ = refine_ms_three_vec;}

  // comparion
  static bool cmpEValueInc(const PrsmPtr &a, const PrsmPtr &b) {
    return a->getEValue() < b->getEValue();
  }

  static bool cmpEValueDec(const PrsmPtr &a, const PrsmPtr &b) {
    return a->getEValue() > b->getEValue();
  }

  static bool cmpMatchFragmentDec(const PrsmPtr &a, const PrsmPtr &b) {
    return a->getMatchFragNum() > b->getMatchFragNum();
  }

  static bool cmpNormMatchFragmentDec(const PrsmPtr &a, const PrsmPtr &b) {
    if (a->getNormMatchFragNum() == b->getNormMatchFragNum()) {
      return a->getProteoformPtr()->getVariablePtmNum() <
          b->getProteoformPtr()->getVariablePtmNum(); 
    } else {
      return a->getNormMatchFragNum() > b->getNormMatchFragNum();
    }
  }

  static bool cmpMatchFragmentDecMatchPeakDec(const PrsmPtr &a, const PrsmPtr &b);

  // sort by number of matched fragment ions, then start position 
  static bool cmpMatchFragDecStartPosInc(const PrsmPtr &a, const PrsmPtr &b);

  // sort by the order of spectrum id, the precursor id
  static bool cmpSpectrumIdIncPrecursorIdInc(const PrsmPtr &a, const PrsmPtr &b);

  // sort by spectrum id then match ions
  static bool cmpSpectrumIdIncMatchFragDec(const PrsmPtr &a, const PrsmPtr &b);

  // sort by spectrum id then evalue
  static bool cmpSpectrumIdIncEvalueInc(const PrsmPtr &a, const PrsmPtr &b);

  // other functions
  xercesc::DOMElement* toXmlElement(XmlDOMDocument* xml_doc);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  void parseXml(xercesc::DOMElement *element);

  static std::string getXmlElementName() {return "prsm";}

  /*
     void initMatchNum(double min_mass);

     bool isMatchMs(MsHeaderPtr header_ptr);
     */
 private:
  int prsm_id_ = -1;
  /* spectrum information */
  int spectrum_id_;

  std::string spectrum_scan_;

  int precursor_id_;

  int spectrum_num_;

  int prec_feature_id_ = -1;

  double prec_feature_inte_ = -1;

  double ori_prec_mass_;
  /* adjusted precursor mass */
  double adjusted_prec_mass_;

  /* protein sequence */
  ProteoformPtr proteoform_ptr_;

  ExtremeValuePtr extreme_value_ptr_;

  double fdr_ = -1;

  double proteoform_fdr_ = -1;

  /* The following are not saved in xml */
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  /* adjusted extended msThree, used for matching ions and peaks */
  ExtendMsPtrVec refine_ms_three_vec_;

  double match_peak_num_ = 0;
  double match_fragment_num_ = 0;

  // functions for initialization
  void init(SpParaPtr sp_para_ptr);
  void initMatchNum(double min_mass);
  void initScores(SpParaPtr sp_para_ptr);
};

typedef std::vector<PrsmPtr> PrsmPtrVec;
typedef std::vector<PrsmPtrVec> PrsmPtrVec2D;
typedef std::vector<PrsmPtrVec2D> PrsmPtrVec3D;

}
#endif

