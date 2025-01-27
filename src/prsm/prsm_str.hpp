//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

  int getProteoClusterId() {return proteo_cluster_id_;}

  int getProtId() {return prot_id_;}

  double getProteoInte() {return proteo_inte_;}

  int getPrecursorId() {return precursor_id_;}

  int getFracFeatureId() {return frac_feature_id_;}

  double getFracFeatureInte() {return frac_feature_inte_;}

  double getFracFeatureApexTime() {return frac_feature_apex_time_;}

  double getFracFeatureMinTime() {return frac_feature_min_time_;}

  double getFracFeatureMaxTime() {return frac_feature_max_time_;}

  int getUnexpectedPtmNum() {return unexpected_ptm_num_;}

  int getVariablePtmNum() {return variable_ptm_num_;}

  std::string getProteoformMatchSeq() {return proteoform_match_seq_;}

  std::string getProteoformDbSeq() {return proteoform_db_seq_;}

  double getMatchPeakNum() {return match_peak_num_;}

  double getMatchFragNum() {return match_frag_num_;}

  double getNormMatchFragNum() {return norm_match_frag_num_;}

  double getEValue() {return e_value_;}

  double getFdr() {return fdr_;}

  double getProteoformFdr() {return proteoform_fdr_;}

  double getOriPrecMass() {return ori_prec_mass_;}

  double getAdjustedPrecMass() {return adjusted_prec_mass_;}

  std::vector<MassShiftPtr> getMassShiftVec() {return mass_shift_vec_;}

  void setFileName(const std::string & fname);

  void setSpectrumId(int id);

  void setProteoClusterId(int id);

  void setProteoInte(double inte);

  void setProtId(int id);

  void setFracFeatureId(int id);

  void setPrecursorId(int id);

  void setFracFeatureInte(double inte);

  void setFracFeatureScore(double score);

  void setFracFeatureApexTime(double apex_time); 

  void setFracFeatureMinTime(double min_time); 

  void setFracFeatureMaxTime(double max_time); 

  void setFdr(double fdr);

  void setProteoformFdr(double proteoform_fdr);

  static bool cmpEValueIncProtInc(const PrsmStrPtr &a, const PrsmStrPtr &b);

  static bool cmpMatchFragDecMatchPeakDecProtInc(const PrsmStrPtr &a, const PrsmStrPtr &b);

  static bool cmpNormMatchFragDecProtInc(const PrsmStrPtr &a, const PrsmStrPtr &b); 

  static bool cmpSpecIncPrecIncEvalueIncProtInc(const PrsmStrPtr &a, const PrsmStrPtr &b);

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

  int precursor_id_;

  int frac_feature_id_;

  double frac_feature_inte_;

  double frac_feature_apex_time_;

  double frac_feature_min_time_;

  double frac_feature_max_time_;

  double ori_prec_mass_;

  double adjusted_prec_mass_;

  // The information from prot_id to mass_shift_vec  
  // is stored in the proteoform class
  int prot_id_;

  std::string seq_name_;

  std::string seq_desc_;

  int proteo_cluster_id_;

  double proteo_inte_;

  int unexpected_ptm_num_;

  int variable_ptm_num_;

  std::string proteoform_match_seq_;

  std::string proteoform_db_seq_;

  int proteoform_start_pos_;

  int proteoform_end_pos_;

  std::vector<MassShiftPtr> mass_shift_vec_;

  //The information below is stored in the prsm class
  double match_peak_num_;

  double match_frag_num_;

  double norm_match_frag_num_;

  double e_value_;

  double fdr_;

  double proteoform_fdr_;


};

typedef std::vector<PrsmStrPtr> PrsmStrPtrVec;
typedef std::vector<PrsmStrPtrVec> PrsmStrPtrVec2D;

}  // namespace toppic

#endif
