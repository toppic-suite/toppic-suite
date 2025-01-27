//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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

#include <cmath>
#include <iomanip>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/base/mod_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "merge/feature_sample_merge.hpp"

namespace toppic {

FeatureSampleMerge::FeatureSampleMerge(const std::vector<std::string> &input_file_names,
                                       const std::string &output_file_name,
                                       const std::string &tool_name, 
                                       double error_tole):
  input_file_names_(input_file_names),
  output_file_name_(output_file_name),
  tool_name_(tool_name),
  error_tole_(error_tole) {}



// normalize by the largest retention time 
void normalizeTime(FeaturePrsmPtrVec &features) {
  double max_time = 1;
  for (size_t i = 0; i < features.size(); i++) {
    if (features[i]->getMaxTime() > max_time) {
      max_time = features[i]->getMaxTime();
    }
  } 
  for (size_t i = 0; i < features.size(); i++) {
    double align_min_time = features[i]->getMinTime()/max_time;
    features[i]->setAlignMinTime(align_min_time);
    double align_max_time = features[i]->getMaxTime()/max_time;
    features[i]->setAlignMaxTime(align_max_time);
    double align_apex_time = features[i]->getApexTime()/max_time;
    features[i]->setAlignApexTime(align_apex_time);
  }
}

class Cell {
 public:
  Cell() {
    sim_ = -1;
    pre_ = -1;
  }
  double sim_ = -1;
  // 0: left, 1: up, 2, diag
  double pre_ = -1;
};

typedef std::shared_ptr<Cell> CellPtr;
typedef std::vector<CellPtr> CellPtrVec;
typedef std::vector<CellPtrVec> CellPtrVec2D;


inline FeaturePrsmPtrVec getTopFeatures(FeaturePrsmPtrVec &features) {
  std::sort(features.begin(), features.end(), FeaturePrsm::cmpInteDec);
  if (features.size() <= 1000) {
    return features;
  }
  FeaturePrsmPtrVec::const_iterator first = features.begin();
  FeaturePrsmPtrVec::const_iterator last = features.begin() + 1000;
  FeaturePrsmPtrVec result(first, last);
  return result;
}

inline CellPtrVec2D initCells(size_t a_size, size_t b_size) {
  CellPtrVec2D cell_2d;
  for (size_t i = 0; i < a_size + 1; i++) {
    CellPtrVec cells(b_size + 1);
    for (size_t j = 0; j < b_size + 1; j++) {
      cells[j] = std::make_shared<Cell>();
    }
    cell_2d.push_back(cells);
  }
  cell_2d[0][0]->sim_ = 0;
  for (size_t i = 1; i < a_size + 1; i++) {
    cell_2d[i][0]->sim_ = 0;
    cell_2d[i][0]->pre_ = 1; // up
  }
  for (size_t j = 1; j < b_size; j++) {
    cell_2d[0][j]->sim_ = 0;
    cell_2d[0][j]->pre_ = 0; // left
  }
  return cell_2d;
}

int findMatch(FeaturePrsmPtr a, FeaturePrsmPtrVec &second, size_t b_index, double mass_tole) {
  double time_tole = 0.1;
  double a_mass = a->getPrecMass();
  //if (mass_tole < 0.01) {
  //  mass_tole = 0.01;
  //}
  double b_time = second[b_index]->getAlignApexTime();
  for (int i = b_index; i >= 0; i--) {
    double mass = second[i]->getPrecMass();
    double time = second[i]->getAlignApexTime();
    if (b_time - time > time_tole) {
      break;
    }
    double error = std::abs(a_mass - mass);
    if (error <= mass_tole) {
      return i;
    }
  }
  for (size_t i = b_index; i < second.size(); i++) {
    double mass = second[i]->getPrecMass();
    double time = second[i]->getAlignApexTime();
    if (time - b_time > time_tole) {
      break;
    }
    double error = std::abs(a_mass - mass);
    if (error <= mass_tole) {
      return i;
    }
  }
  return -1;
}


void sampleAlignTime(FeaturePrsmPtrVec &first, FeaturePrsmPtrVec &second, 
                     double mass_tole, std::vector<std::pair<double,double>> &time_pairs) {
  FeaturePrsmPtrVec a = getTopFeatures(first);
  std::sort(a.begin(), a.end(), FeaturePrsm::cmpTimeInc);
  FeaturePrsmPtrVec b = getTopFeatures(second);
  std::sort(b.begin(), b.end(), FeaturePrsm::cmpTimeInc);

  size_t a_size = a.size();
  size_t b_size = b.size();
  CellPtrVec2D cell_2d = initCells(a_size, b_size);
  //alignment
  for (size_t i = 1; i < a_size + 1; i++) {
    for (size_t j = 1; j < b_size + 1; j++) {
      double left_sim = cell_2d[i][j-1]->sim_;
      double up_sim = cell_2d[i-1][j]->sim_;
      double diag_sim = cell_2d[i-1][j-1]->sim_;
      int pos = findMatch(a[i-1], b, j-1, mass_tole); 
      // if there is a match
      if (pos >= 0) {
        double inte = a[i-1]->getProteoInte();
        double log_inte =  std::log10(inte + 1);
        diag_sim += log_inte;
      }
      if (diag_sim >= left_sim && diag_sim >= up_sim) {
        cell_2d[i][j]->sim_ = diag_sim;
        cell_2d[i][j]->pre_ = 2; // diag
      }
      else if (up_sim >= left_sim) {
        cell_2d[i][j]->sim_ = up_sim;
        cell_2d[i][j]->pre_ = 1; // up
      }
      else {
        cell_2d[i][j]->sim_ = left_sim;
        cell_2d[i][j]->pre_ = 0; // left
      }
    }
  }
  /*
  for (size_t i = 0; i < a_size; i++) {
    LOG_ERROR("A " << i << " time " << a[i]->getAlignApexTime());
  }
  for (size_t i = 0; i < b_size; i++) {
    LOG_ERROR("B " << i << " time " << b[i]->getAlignApexTime());
  }
  */
  //backtracking
  size_t i = a_size;
  size_t j = b_size;
  std::pair<double,double> last_pair(1,1);
  time_pairs.push_back(last_pair);
  while (i > 0 && j > 0) {
    int pre = cell_2d[i][j]->pre_;
    if (pre == 2) {// diag
      int pos = findMatch(a[i-1], b, j-1, mass_tole); 
      if (pos >= 0) {
        double a_time = a[i-1]->getAlignApexTime();
        double b_time = b[pos]->getAlignApexTime();
        std::pair<double, double>time_pair(a_time,b_time);
        time_pairs.push_back(time_pair);
        LOG_DEBUG("i " << (i-1) 
                  << " time " << a[i-1]->getAlignApexTime() 
                  << " mass " << a[i-1]->getPrecMass() 
                  << " intensity " << a[i-1]->getProteoInte() 
                  << " j " << (j-1) 
                  << " j " << pos 
                  << " time " << b[pos]->getAlignApexTime() 
                  << " mass " << b[pos]->getPrecMass() 
                  << " intensity " << b[pos]->getProteoInte()); 
      }
      i--;
      j--;
    }
    else if (pre == 1) { // up
      i--;
    }
    else {
      j--;
    }
  }
  std::pair<double,double> first_pair(0,0);
  time_pairs.push_back(first_pair);
  std::reverse(time_pairs.begin(), time_pairs.end());
}

double findAlignTime(double time, std::vector<std::pair<double,double>> &time_pairs) {
  for (size_t i = 0; i < time_pairs.size() - 1; i++) {
    double a = time_pairs[i].second;
    double b = time_pairs[i+1].second;
    LOG_DEBUG("a " << a  << " b " << b << " time " << time);
    if (time >= a && time <= b) {
      double align_a = time_pairs[i].first;
      double align_b = time_pairs[i+1].first;
      if (a == b) {
        return align_a;
      }
      else {
        double align_time  = align_a + (align_b - align_a) * (time-a)/(b-a);
        LOG_DEBUG("time " << time   << " align time  " << align_time );
        return align_time;
      }
    }
  }
  LOG_ERROR("Can not find align time!");
  return -1;
}

void setAlignTime(FeaturePrsmPtrVec &features, 
                  std::vector<std::pair<double,double>> &time_pairs) {
  for (size_t i = 0; i < features.size(); i++) {
    double time_begin = findAlignTime(features[i]->getAlignMinTime(), time_pairs);
    features[i]->setAlignMinTime(time_begin);
    double time_end = findAlignTime(features[i]->getAlignMaxTime(), time_pairs);
    features[i]->setAlignMaxTime(time_end);
    double time_apex = findAlignTime(features[i]->getAlignApexTime(), time_pairs);
    LOG_DEBUG("ori align time " << features[i]->getAlignApexTime() << " normalized " << time_apex); 
    features[i]->setAlignApexTime(time_apex);
  }
}

//sort feature 2d by intensities
FeaturePrsmPtr findMatchFeature(FeaturePrsmPtr feature, FeaturePrsmPtrVec &features, 
                                double mass_tole) {
  double time_tole = 0.3;
  double strict_mass_tole = 0.01;
  double strict_time_tole = 0.01;
  double cur_mass = feature->getPrecMass();
  double cur_time = feature->getAlignApexTime();
  std::string cur_prot = feature->getProtName();
  for (size_t i = 0; i < features.size(); i++) {
    if (features[i] == nullptr) {
      continue;
    }
    std::string prot = features[i]->getProtName();
    double mass = features[i]->getPrecMass();
    double time = features[i]->getAlignApexTime();
    if (cur_prot == prot) {
      if (abs(time-cur_time) <= time_tole && abs(mass - cur_mass) <= mass_tole) {
        FeaturePrsmPtr result = features[i];
        features[i] = nullptr;
        return result;
      }
    }
    if (prot == "") {
      if (abs(time-cur_time) <= strict_time_tole 
          && abs(mass - cur_mass) <= strict_mass_tole) {
        FeaturePrsmPtr result = features[i];
        features[i] = nullptr;
        return result;
      }
    }

  }
  return nullptr;
}

bool featureUsed(FeaturePrsmPtr feature, FeaturePrsmPtrVec2D &table) {
  int cur_sample = feature->getSampleId();
  int feature_id = feature->getProteoId();
  for (size_t i = 0; i < table.size(); i++) {
    FeaturePrsmPtr cur_feature= table[i][cur_sample];
    if (cur_feature != nullptr) {
      int cur_id = cur_feature->getProteoId();
      if (feature_id == cur_id) {
        return true;
      }
    }
  }
  return false;
}

void getFeatureTable(FeaturePrsmPtrVec& all_features, FeaturePrsmPtrVec2D& features_2d, 
                     FeaturePrsmPtrVec2D& table, FeaturePrsmPtrVec &examples, double mass_tole) {
  size_t sample_size = features_2d.size();
  for (size_t i = 0; i < all_features.size(); i++) {
    if (all_features[i]->getProtName() == "") {
      continue;
    }
    // to do check if it is used
    if (featureUsed(all_features[i], table)) {
      continue;
    }
    int cur_sample = all_features[i]->getSampleId();

    FeaturePrsmPtrVec cur_row(sample_size);
    for (size_t j = 0; j < sample_size; j++) {
      if ((int)j == cur_sample) {
        cur_row[j] = all_features[i];
      }
      else {
        FeaturePrsmPtr feature_ptr = findMatchFeature(all_features[i], features_2d[j], mass_tole);
        if (feature_ptr != nullptr && !featureUsed(feature_ptr, table)) {
          cur_row[j] = feature_ptr;
        }
      }
    }
    table.push_back(cur_row);
    examples.push_back(all_features[i]);
  }
}

void FeatureSampleMerge::outputTable(FeaturePrsmPtrVec2D &table, 
                                     FeaturePrsmPtrVec &examples, int sample_num) {
  std::ofstream file;
  file.open(output_file_name_.c_str());
  // write title
  std::string delim = "\t";
  file << "Protein accession" << delim
    << "Protein description" << delim
    << "First residue" << delim
    << "Last residue" << delim
    << "Proteoform" << delim
    << "Precursor mass" << delim
    << "# matched samples" << delim;
    //<< "# matched samples with MS/MS ID" << delim;

  for (int i = 0; i < sample_num; i++) {
    file 
      << input_file_names_[i] << " Abundance" << delim
      << input_file_names_[i] << " Spectrum id" << delim
      << input_file_names_[i] << " Retention time begin" << delim
      << input_file_names_[i] << " Retention time end" << delim
      << input_file_names_[i] << " Normalized time apex" << delim;
  }
  file << std::endl;
  int cluster_num = table.size();
  int match_num = 0;
  for (int i = 0; i < cluster_num; i++) {
    FeaturePrsmPtr feature_ptr = examples[i];
    file << feature_ptr->getProtName() << delim
      << "\"" << feature_ptr->getProtDesc() << "\"" << delim
      << (feature_ptr->getFirstResidue() + 1) << delim
      << (feature_ptr->getLastResidue() + 1) << delim
      << feature_ptr->getProteoform() << delim
      << feature_ptr->getPrecMass() << delim;
    int sample_match_num = 0; 
    int sample_id_match_num = 0;
    for (int j = 0; j < sample_num; j++) {
      if (table[i][j] != nullptr) {
        sample_match_num ++;
        if (table[i][j]->getMs2Id() >= 0) {
          sample_id_match_num ++;
        }
      }
    }
    file << sample_match_num << delim;
    LOG_DEBUG("Sample ID match number: " << sample_id_match_num); 
    if (sample_match_num == sample_num) {
      match_num = match_num + 1;
    }
    for (int j = 0; j < sample_num; j++) {
      FeaturePrsmPtr sample_feature = table[i][j];
      if (sample_feature == nullptr) {
        file << delim << delim << delim << delim << delim;
      }
      else {
        /*
        if (sample_feature->getMs2Id()>= 0) {
          file << "Match with MS/MS ID" << delim;
        }
        else {
          file << "Match without MS/MS ID" << delim; 
        }
        */
        file <<  std::setprecision(3) << std::scientific << sample_feature->getProteoInte() << delim;
        if (sample_feature->getMs2Id()>= 0) {
          file << std::fixed << sample_feature->getMs2Id() << delim;
        }
        else {
          file << std::fixed << delim;
        }
        file << sample_feature->getMinTime() << delim
          << sample_feature->getMaxTime() << delim
          << sample_feature->getAlignApexTime() << delim;
      }
    }
    file << std::endl;
  }
  file.close();
  std::cout << "Proteoform number:" << cluster_num << " matched number: " << match_num << std::endl;
}


void FeatureSampleMerge::process() {
  FeaturePrsmPtrVec2D features_2d; 
  FeaturePrsmPtrVec all_features; 
  size_t sample_num = input_file_names_.size();

  for (size_t k = 0; k < sample_num; k++) {
    std::string input_file_name = input_file_names_[k];
    std::string base_name = file_util::basename(input_file_name);
    if (str_util::endsWith(base_name, "_ms2")) {
      base_name = base_name.substr(0, base_name.size() - 4);
    }
    else {
      LOG_ERROR("The file name " << input_file_name << " does not end with _ms2.msalign!");
    }

    std::string prsm_file_name = base_name + "_ms2_" + tool_name_ + "_proteoform.xml";
    PrsmStrPtrVec prsms = prsm_reader_util::readAllPrsmStrsMatchSeq(prsm_file_name);
    if (prsms.size() == 0) {
      LOG_WARN("The file " << prsm_file_name  << " does not contain any PrSM identifications!");
    }
    FeaturePrsmPtrVec features;
    for (size_t i = 0; i < prsms.size(); i++) {
      FeaturePrsmPtr feat_ptr = std::make_shared<FeaturePrsm>(prsms[i]);
      features.push_back(feat_ptr);
    }
    std::sort(features.begin(), features.end(), FeaturePrsm::cmpTimeInc);
    for (size_t i = 0; i < features.size(); i++) {
      features[i]->setSampleId(k);
    }
    normalizeTime(features);
    if (k > 0) {
      std::vector<std::pair<double,double>> time_pairs;
      sampleAlignTime(features_2d[0], features, error_tole_, time_pairs);
      setAlignTime(features, time_pairs);
    }
    features_2d.push_back(features);
    all_features.insert(all_features.end(), features.begin(), features.end());
  }
  std::sort(all_features.begin(), all_features.end(), FeaturePrsm::cmpInteDec);
  for (size_t k = 0; k < sample_num; k++) {
    std::sort(features_2d[k].begin(), features_2d[k].end(), FeaturePrsm::cmpInteDec);
  }
  FeaturePrsmPtrVec2D table; 
  FeaturePrsmPtrVec examples; 
  getFeatureTable(all_features, features_2d, table, examples, error_tole_);
  outputTable(table, examples, sample_num);
}

}  // namespace prot


