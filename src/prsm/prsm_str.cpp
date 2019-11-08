//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#include <limits>
#include <string>
#include <vector>
#include <algorithm>

#include "common/util/logger.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_str.hpp"

namespace toppic {

PrsmStr::PrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  std::string line = prsm_util::getXmlLine(str_vec_, "<spectrum_id>");
  file_name_ = prsm_util::getValueStr(line);
  line = prsm_util::getXmlLine(str_vec_, "<spectrum_id>");
  spectrum_id_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<spectrum_scan>");
  spectrum_scan_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<precursor_id>");
  precursor_id_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<precursor_feature_id>");
  precursor_feature_id_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<precursor_feature_inte>");
  precursor_feature_inte_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<ori_prec_mass>");
  ori_prec_mass_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<adjusted_prec_mass>");
  adjusted_prec_mass_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<seq_name>");
  seq_name_ = prsm_util::getValueStr(line);
  line = prsm_util::getXmlLine(str_vec_, "<seq_desc>");
  seq_desc_ = prsm_util::getValueStr(line);
  line = prsm_util::getXmlLine(str_vec_, "<match_fragment_num>");
  match_frag_num_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<norm_match_fragment_num>");
  norm_match_frag_num_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<e_value>");
  if (line == "") {
    e_value_ = 0.0;
  } else {
    std::string str = prsm_util::getValueStr(line);
    LOG_DEBUG("e value string " << str);
    e_value_ = str_util::scientificToDouble(str);
    LOG_DEBUG("e value value " << e_value_);
  }
  line = prsm_util::getXmlLine(str_vec_, "<fdr>");
  fdr_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<proteoform_fdr>");
  proteoform_fdr_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<start_pos>");
  proteoform_start_pos_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<end_pos>");
  proteoform_end_pos_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<proteo_cluster_id>");
  cluster_id_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<prot_id>");
  prot_id_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<unexpected_ptm_num>");
  unexpected_ptm_num_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<variable_ptm_num>");
  variable_ptm_num_ = std::stoi(prsm_util::getValueStr(line));

  std::vector<std::string> mass_lines = prsm_util::getXmlLineVec(str_vec_, "<shift>");
  std::vector<std::string> left_pos_lines = prsm_util::getXmlLineVec(str_vec_, "<left_bp_pos>");
  std::vector<std::string> right_pos_lines = prsm_util::getXmlLineVec(str_vec_, "<right_bp_pos>");

  for (size_t i = 0; i < mass_lines.size(); i++) {
    mass_shift_vec_.push_back(std::make_shared<MassShiftStr>(std::stod(prsm_util::getValueStr(mass_lines[i])),
                                                             std::stoi(prsm_util::getValueStr(left_pos_lines[i])),
                                                             std::stoi(prsm_util::getValueStr(right_pos_lines[i]))));
  }
}

int getXmlLineIndex(const std::vector<std::string> &str_vec,
                    const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return i;
    }
  }
  return -1;
}

bool PrsmStr::cmpSpectrumIdIncPrecursorIdInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
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

void PrsmStr::setFdr(double fdr) {
  int i = getXmlLineIndex(str_vec_, "fdr");
  str_vec_[i] = "<fdr>" + str_util::toString(fdr) + "</fdr>";
  fdr_ = fdr;
}

void PrsmStr::setProteoformFdr(double proteoform_fdr) {
  int i = getXmlLineIndex(str_vec_, "proteoform_fdr");
  str_vec_[i] = "<proteoform_fdr>" + str_util::toString(proteoform_fdr) + "</proteoform_fdr>";
  proteoform_fdr_ = proteoform_fdr;
}

void PrsmStr::setFileName(const std::string & fname) {
  int i = getXmlLineIndex(str_vec_, "file_name");
  str_vec_[i] = "<file_name>" + fname + "</file_name>";
  file_name_ = fname;
}

void PrsmStr::setSpectrumId(int id) {
  int i = getXmlLineIndex(str_vec_, "spectrum_id");
  str_vec_[i] = "<spectrum_id>" + str_util::toString(id) + "</spectrum_id>";
  spectrum_id_ = id;
}

void PrsmStr::setPrecFeatureId(int id) {
  int i = getXmlLineIndex(str_vec_, "precursor_feature_id");
  str_vec_[i] = "<precursor_feature_id>" + str_util::toString(id) + "</precursor_feature_id>";
  precursor_feature_id_ = id;
}

void PrsmStr::setPrecFeatureInte(double inte) {
  int i = getXmlLineIndex(str_vec_, "precursor_feature_inte");
  str_vec_[i] = "<precursor_feature_inte>" + str_util::toString(inte) + "</precursor_feature_inte>";
  precursor_feature_inte_ = inte;
}

void PrsmStr::setFracFeatureScore(double score) {
  int i = getXmlLineIndex(str_vec_, "frac_feature_score");
  str_vec_[i] = "<frac_feature_score>" + str_util::toString(score) + "</frac_feature_score>";
}


void PrsmStr::setPrecursorId(int id) {
  int i = getXmlLineIndex(str_vec_, "precursor_id");
  str_vec_[i] = "<precursor_id>" + str_util::toString(id) + "</precursor_id>";
  precursor_id_ = id;
}

void PrsmStr::setClusterId(int id) {
  int i = getXmlLineIndex(str_vec_, "proteo_cluster_id");
  str_vec_[i] = "<proteo_cluster_id>" + str_util::toString(id) + "</proteo_cluster_id>";
  cluster_id_ = id;
}

void PrsmStr::setProtId(int id) {
  int i = getXmlLineIndex(str_vec_, "prot_id");
  str_vec_[i] = "<prot_id>" + str_util::toString(id) + "</prot_id>";
  prot_id_ = id;
}

bool PrsmStr::isSameSeqAndMass(const PrsmStrPtr &a, const PrsmStrPtr &b, double ppo) {
  if (a->getSeqName() != b->getSeqName()) {
    return false;
  }
  if (a->getProteoformStartPos() != b->getProteoformStartPos()) {
    return false;
  }
  if (a->getProteoformEndPos() != b->getProteoformEndPos()) {
    return false;
  }
  double thresh = a->getAdjustedPrecMass() * ppo;
  if (std::abs(a->getAdjustedPrecMass() - b->getAdjustedPrecMass()) > thresh) {
    return false;
  }
  return true;
}

bool PrsmStr::isStrictCompatiablePtmSpecies(const PrsmStrPtr & a, const PrsmStrPtr & b, double ppo) {
  if (!isSameSeqAndMass(a, b, ppo)) {
    return false;
  }

  if (a->getChangeStrVec().size() != b->getChangeStrVec().size()) {
    return false;
  }

  double shift_tolerance = a->getAdjustedPrecMass() * ppo;
  std::vector<MassShiftStrPtr> a_shift_vec = a->getChangeStrVec();
  std::vector<MassShiftStrPtr> b_shift_vec = b->getChangeStrVec();
  std::sort(a_shift_vec.begin(), a_shift_vec.end(), MassShiftStr::cmpPosInc);
  std::sort(b_shift_vec.begin(), b_shift_vec.end(), MassShiftStr::cmpPosInc);
  for (size_t i = 0; i < a->getChangeStrVec().size(); i++) {
    MassShiftStrPtr ac = a_shift_vec[i];
    MassShiftStrPtr bc = b_shift_vec[i];
    if (ac->right_pos_ <= bc->left_pos_ || bc->right_pos_ <= ac->left_pos_) {
      return false;
    }
    if (std::abs(ac->mass_shift_ - bc->mass_shift_) > shift_tolerance) {
      return false;
    }
  }
  return true;
}

}  // namespace toppic
