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

#ifndef TOPPIC_PRSM_PRSM_STR_HPP_
#define TOPPIC_PRSM_PRSM_STR_HPP_

#include <vector>
#include <string>

#include "seq/mass_shift.hpp"

namespace toppic {

class PrsmStr;
typedef std::shared_ptr<PrsmStr> PrsmStrPtr;

class PrsmStr {
 public:
  explicit PrsmStr(const std::vector<std::string> &str_vec);

  std::vector<std::string> getStrVec() {return str_vec_;}

  std::string getFileName() {return file_name_;}

  int getSpectrumId() {return spectrum_id_;}
  
  int getSpectrumScan() {return spectrum_scan_;}

  std::string getSeqName() {return seq_name_;}

  std::string getSeqDesc() {return seq_desc_;}

  int getProteoformStartPos() {return proteoform_start_pos_;}

  int getProteoformEndPos() {return proteoform_end_pos_;}

  int getClusterId() {return cluster_id_;}

  int getProtId() {return prot_id_;}

  int getPrecursorId() {return precursor_id_;}

  int getSampleFeatureId() {return sample_feature_id_;}

  int getUnexpectedPtmNum() {return unexpected_ptm_num_;}

  int getVariablePtmNum() {return variable_ptm_num_;}

  std::string getProteoformMatchSeq() {return proteoform_match_seq_;}

  double getPrecFeatureInte() {return sample_feature_inte_;}

  double getMatchFragNum() {return match_frag_num_;}

  double getNormMatchFragNum() {return norm_match_frag_num_;}

  double getEValue() {return e_value_;}

  double getFdr() {return fdr_;}

  double getProteoformFdr() {return proteoform_fdr_;}

  double getOriPrecMass() {return ori_prec_mass_;}

  double getAdjustedPrecMass() {return adjusted_prec_mass_;}

  std::vector<MassShiftPtr> getMassShiftVec() {return mass_shift_vec_;}

  //std::string getProteinMatchSeq() {return protein_match_seq_;}

  int getSampleId() {return sample_id_;}

  double getTimeApex() {return time_apex_;}

  void setFileName(const std::string & fname);

  void setSpectrumId(int id);

  void setClusterId(int id);

  void setProtId(int id);

  void setPrecFeatureId(int id);

  void setPrecursorId(int id);

  void setPrecFeatureInte(double inte);

  void setFracFeatureScore(double score);

  void setFdr(double fdr);

  void setProteoformFdr(double proteoform_fdr);

  //void setProteinMatchSeq(const std::string & seq) {protein_match_seq_ = seq;}

  void setSampleId(int sample_id) {sample_id_ = sample_id;}

  void setTimeApex(double time_apex) {time_apex_ = time_apex;}

  static bool cmpEValueInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getEValue() < b->getEValue();}

  static bool cmpMatchFragmentDec(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getMatchFragNum() > b->getMatchFragNum();}

  static bool cmpNormMatchFragmentDec(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getNormMatchFragNum() > b->getNormMatchFragNum();}

  static bool cmpSpectrumIdInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getSpectrumId() < b->getSpectrumId();}

  static bool cmpSpectrumIdIncPrecursorIdInc(const PrsmStrPtr &a, const PrsmStrPtr &b);

  static bool isSameSeq(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getSeqName() == b->getSeqName();
  }

  static bool isSameSeqAndMass(const PrsmStrPtr &a, const PrsmStrPtr &b, double ppo);

  static bool isSimpleMatch(const PrsmStrPtr &a, const PrsmStrPtr &b, double ppo);

  static bool isStrictCompatiablePtmSpecies(const PrsmStrPtr & a, const PrsmStrPtr & b, double ppo);

 private:
  std::vector<std::string> str_vec_;

  std::string file_name_;

  int spectrum_id_;

  int spectrum_scan_;

  std::string seq_name_;

  std::string seq_desc_;

  int cluster_id_;

  int prot_id_;

  int precursor_id_;

  int sample_feature_id_;

  double sample_feature_inte_;

  int unexpected_ptm_num_;

  int variable_ptm_num_;

  std::string proteoform_match_seq_;

  int proteoform_start_pos_;

  int proteoform_end_pos_;

  double ori_prec_mass_;

  double adjusted_prec_mass_;

  double match_frag_num_;

  double norm_match_frag_num_;

  double e_value_;

  double fdr_;

  double proteoform_fdr_;

  std::vector<MassShiftPtr> mass_shift_vec_;

  //std::string protein_match_seq_;

  int sample_id_;

  double time_apex_ = -1;
};

typedef std::vector<PrsmStrPtr> PrsmStrPtrVec;
typedef std::vector<PrsmStrPtrVec> PrsmStrPtrVec2D;

}  // namespace toppic

#endif
