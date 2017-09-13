// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_PRSM_PRSM_STR_HPP_
#define PROT_PRSM_PRSM_STR_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class PrsmStr;
typedef std::shared_ptr<PrsmStr> PrsmStrPtr;

class PrsmStr {
 public:
  explicit PrsmStr(const std::vector<std::string> &str_vec);

  std::vector<std::string> getStrVec() {return str_vec_;}

  int getSpectrumId() {return spectrum_id_;}

  std::string getSeqName() {return seq_name_;}

  std::string getSeqDesc() {return seq_desc_;}

  int getSpeciesId() {return species_id_;}

  int getProtId() {return prot_id_;}

  int getPrecursorId() {return precursor_id_;}

  int getPrecFeatureId() {return precursor_feature_id_;}

  double getPrecFeatureInte() {return precursor_feature_inte_;}

  double getMatchFragNum() {return match_frag_num_;}

  double getNormMatchFragNum() {return norm_match_frag_num_;}

  double getEValue() {return e_value_;}

  double getFdr() {return fdr_;}

  double getProteoformFdr() {return proteoform_fdr_;}

  double getAdjustedPrecMass() {return adjusted_prec_mass_;}

  void setSpectrumId(int id);

  void setSpeciesId(int id);

  void setProtId(int id);

  void setPrecFeatureId(int id);

  void setPrecFeatureInte(double inte);

  void setFdr(double fdr);

  void setProteoformFdr(double proteoform_fdr);

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

 private:
  std::vector<std::string> str_vec_;

  int spectrum_id_;

  std::string seq_name_;

  std::string seq_desc_;

  int species_id_;

  int prot_id_;

  int precursor_id_;

  int precursor_feature_id_;

  double precursor_feature_inte_;

  double adjusted_prec_mass_;

  double match_frag_num_;

  /* a tempory variable for testing mass graph alignment */
  double norm_match_frag_num_;

  double e_value_;

  double fdr_;

  double proteoform_fdr_;
};

typedef std::vector<PrsmStrPtr> PrsmStrPtrVec;
typedef std::vector<PrsmStrPtrVec> PrsmStrPtr2D;

}  // namespace prot

#endif
