//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_PRSM_PRSM_HPP_
#define TOPPIC_PRSM_PRSM_HPP_

#include <string>
#include <vector>

#include "seq/proteoform.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/extend_ms.hpp"
#include "para/sp_para.hpp"
#include "prsm/expected_value.hpp"

namespace toppic {

class Prsm;
typedef std::shared_ptr<Prsm> PrsmPtr;

class Prsm {
 public:
  Prsm(ProteoformPtr proteoform_ptr, const DeconvMsPtrVec &deconv_ms_ptr_vec,
       double adjusted_prec_mass, SpParaPtr sp_para_ptr);

  Prsm(XmlDOMElement* element, FastaIndexReaderPtr reader_ptr,
       const ModPtrVec &fix_mod_list);

  std::string getFileName() {return file_name_;}

  int getPrsmId() {return prsm_id_;}

  int getSpectrumId() {return spectrum_id_;}

  std::string getSpectrumScan() {return spectrum_scan_;}

  int getPrecursorId() {return precursor_id_;}

  int getSampleFeatureId() {return sample_feature_id_;}

  double getSampleFeatureInte() {return sample_feature_inte_;}

  double getFracFeatureScore() {return frac_feature_score_;}

  double getOriPrecMass() {return ori_prec_mass_;}

  double getAdjustedPrecMass() {return adjusted_prec_mass_;}

  ProteoformPtr getProteoformPtr() {return proteoform_ptr_;}

  ExpectedValuePtr getExpectedValuePtr() {return expected_value_ptr_;}

  double getFdr() {return fdr_;}

  double getProteoformFdr() {return proteoform_fdr_;}

  DeconvMsPtrVec getDeconvMsPtrVec() {return deconv_ms_ptr_vec_;}

  ExtendMsPtrVec getRefineMsPtrVec() {return refine_ms_three_vec_;}

  double getMatchPeakNum() {return match_peak_num_;}

  double getMatchFragNum() {return match_fragment_num_;}

  // Normalized match fragment number is used in TopMG for ranking PrSMs
  double getNormMatchFragNum();

  // Expected related functions
  double getEValue();

  double getPValue();

  double getOneProtProb();

  double getTimeApex() {return time_apex_;};

  int getHitCnt() {return hit_cnt_;}

  bool getIsExactMatch() {return is_exact_match_;}

  XmlDOMElement* getElement() {return element_;}

  FastaIndexReaderPtr getReaderPtr() {return reader_ptr_;}

  const ModPtrVec& getFixModList() {return fix_mod_list_;}

  // setter
  void setFileName(const std::string & fname) {file_name_ = fname;}

  void setPrsmId(int id) {prsm_id_ = id;}

  void setSpectrumId(int spectrum_id) {spectrum_id_ = spectrum_id;}

  void setSpectrumScan(std::string spectrum_scan) {spectrum_scan_ = spectrum_scan;}

  void setPrecurorId(int precursor_id) {precursor_id_ = precursor_id;}

  void setFracFeatureScore(double score) {frac_feature_score_ = score;}

  void setOriPrecMass(double prec_mass) {ori_prec_mass_ = prec_mass;}

  void setProteoformPtr(ProteoformPtr proteoform) {proteoform_ptr_ = proteoform;}

  void setProteoformPtr(ProteoformPtr proteoform, SpParaPtr sp_para_ptr);

  void setExpectedValuePtr(ExpectedValuePtr ev_ptr) {expected_value_ptr_ = ev_ptr;}

  void setFdr(double fdr) {fdr_ = fdr;}

  void setProteoformFdr(double proteoform_fdr) {proteoform_fdr_ = proteoform_fdr;}

  void setDeconvMsPtrVec(DeconvMsPtrVec ms_vec) {deconv_ms_ptr_vec_ = ms_vec;}

  void setRefineMsVec(ExtendMsPtrVec refine_ms_three_vec) {
    refine_ms_three_vec_ = refine_ms_three_vec;}

  void setAdjustedPrecMass(double new_prec_mass);

  void setHitCnt(int hit_cnt) {hit_cnt_ = hit_cnt;}

  void setIsExactMatch(bool is_exact_match) {is_exact_match_ = is_exact_match;}

  // comparion
  static bool cmpEValueInc(const PrsmPtr &a, const PrsmPtr &b) {
    return a->getEValue() < b->getEValue();}

  static bool cmpEValueDec(const PrsmPtr &a, const PrsmPtr &b) {
    return a->getEValue() > b->getEValue();}

  static bool cmpMatchFragmentDec(const PrsmPtr &a, const PrsmPtr &b) {
    return a->getMatchFragNum() > b->getMatchFragNum();}

  static bool cmpNormMatchFragmentDec(const PrsmPtr &a, const PrsmPtr &b);

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
  XmlDOMElement* toXmlElement(XmlDOMDocument* xml_doc);

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  void parseXml(XmlDOMElement *element);

  static std::string getXmlElementName() {return "prsm";}

 private:
  std::string file_name_;

  int prsm_id_ = -1;
  /* spectrum information */
  int spectrum_id_;

  std::string spectrum_scan_;

  int precursor_id_;

  int spectrum_num_;

  int sample_feature_id_ = -1;

  double sample_feature_inte_ = -1;

  double frac_feature_score_ = -1000;

  double ori_prec_mass_;
  /* adjusted precursor mass */
  double adjusted_prec_mass_;

  /* protein sequence */
  ProteoformPtr proteoform_ptr_;

  ExpectedValuePtr expected_value_ptr_;

  double fdr_ = -1;

  double proteoform_fdr_ = -1;

  double time_apex_ = -1;

  int hit_cnt_ = 0;

  bool is_exact_match_ = true; 

  /* The following are not saved in xml */
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  /* adjusted extended msThree, used for matching ions and peaks */
  ExtendMsPtrVec refine_ms_three_vec_;

  double match_peak_num_ = 0;
  double match_fragment_num_ = 0;

  XmlDOMElement* element_;
  FastaIndexReaderPtr reader_ptr_;
  ModPtrVec fix_mod_list_;

  // functions for initialization
  void init(SpParaPtr sp_para_ptr);
  void initMatchNum(double min_mass);
  void initScores(SpParaPtr sp_para_ptr);
};

typedef std::vector<PrsmPtr> PrsmPtrVec;
typedef std::vector<PrsmPtrVec> PrsmPtrVec2D;
typedef std::vector<PrsmPtrVec2D> PrsmPtrVec3D;

}  // namespace toppic
#endif

