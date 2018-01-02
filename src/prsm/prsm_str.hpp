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


#ifndef PROT_PRSM_PRSM_STR_HPP_
#define PROT_PRSM_PRSM_STR_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class ChangeStr {
 public:
  ChangeStr(double mass_shift, int left_pos, int right_pos):
      mass_shift_(mass_shift),
      left_pos_(left_pos),
      right_pos_(right_pos) {}

  double mass_shift_;
  int left_pos_, right_pos_;
  static bool cmpPosInc(const std::shared_ptr<ChangeStr> &a,
                        const std::shared_ptr<ChangeStr> &b);
};

typedef std::shared_ptr<ChangeStr> ChangeStrPtr;

class PrsmStr;
typedef std::shared_ptr<PrsmStr> PrsmStrPtr;

class PrsmStr {
 public:
  explicit PrsmStr(const std::vector<std::string> &str_vec);

  std::vector<std::string> getStrVec() {return str_vec_;}

  std::string getFileName() {return file_name_;}

  int getSpectrumId() {return spectrum_id_;}

  std::string getSeqName() {return seq_name_;}

  std::string getSeqDesc() {return seq_desc_;}

  int getProteoformStartPos() {return proteoform_start_pos_;}

  int getProteoformEndPos() {return proteoform_end_pos_;}

  int getClusterId() {return cluster_id_;}

  int getProtId() {return prot_id_;}

  int getPrecursorId() {return precursor_id_;}

  int getPrecFeatureId() {return precursor_feature_id_;}

  int getUnexpectedPtmNum() {return unexpected_ptm_num_;}

  double getPrecFeatureInte() {return precursor_feature_inte_;}

  double getMatchFragNum() {return match_frag_num_;}

  double getNormMatchFragNum() {return norm_match_frag_num_;}

  double getEValue() {return e_value_;}

  double getFdr() {return fdr_;}

  double getProteoformFdr() {return proteoform_fdr_;}

  double getAdjustedPrecMass() {return adjusted_prec_mass_;}

  std::vector<std::shared_ptr<ChangeStr> > getChangeStrVec() {return change_vec_;}

  std::string getProteinMatchSeq() {return protein_match_seq_;}

  void setFileName(const std::string & fname);

  void setSpectrumId(int id);

  void setClusterId(int id);

  void setProtId(int id);

  void setPrecFeatureId(int id);

  void setPrecursorId(int id);

  void setPrecFeatureInte(double inte);

  void setFdr(double fdr);

  void setProteoformFdr(double proteoform_fdr);

  void setProteinMatchSeq(const std::string & seq) {protein_match_seq_ = seq;}

  static bool cmpEValueInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getEValue() < b->getEValue();
  }

  static bool cmpMatchFragmentDec(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getMatchFragNum() > b->getMatchFragNum();
  }

  static bool cmpNormMatchFragmentDec(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getNormMatchFragNum() > b->getNormMatchFragNum();
  }

  static bool cmpSpectrumIdInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getSpectrumId() < b->getSpectrumId();
  }

  static bool cmpSpectrumIdIncPrecursorIdInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    if (a->getSpectrumId() < b->getSpectrumId()) {
      return true;
    } else if (a->getSpectrumId() > b->getSpectrumId()) {
      return false;
    } else {
      if (a->getPrecursorId() < b->getPrecursorId()) {
        return true;
      }
      return false;
    }
  }

  static bool isSameSeqAndMass(const PrsmStrPtr &a, const PrsmStrPtr &b, double ppo);

  static bool isStrictCompatiablePtmSpecies(const PrsmStrPtr & a, const PrsmStrPtr & b, double ppo);

 private:
  std::vector<std::string> str_vec_;

  std::string file_name_;

  int spectrum_id_;

  std::string seq_name_;

  std::string seq_desc_;

  int cluster_id_;

  int prot_id_;

  int precursor_id_;

  int precursor_feature_id_;

  int unexpected_ptm_num_;

  int proteoform_start_pos_;

  int proteoform_end_pos_;

  double precursor_feature_inte_;

  double adjusted_prec_mass_;

  double match_frag_num_;

  double norm_match_frag_num_;

  double e_value_;

  double fdr_;

  double proteoform_fdr_;

  std::vector<std::shared_ptr<ChangeStr> > change_vec_;

  std::string protein_match_seq_;
};

typedef std::vector<PrsmStrPtr> PrsmStrPtrVec;
typedef std::vector<PrsmStrPtrVec> PrsmStrPtr2D;

}  // namespace prot

#endif
