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

#include <iomanip>
#include <algorithm>

#include "common/util/logger.hpp"
#include "deconv/feature/feature_prsm_reader.hpp"
#include "deconv/feature/feature_sample_merge.hpp"

namespace toppic {

FeatureSampleMerge::FeatureSampleMerge(const std::vector<std::string> &input_file_names,
                                       const std::string &output_file_name,
                                       double ppm):
    input_file_names_(input_file_names),
    output_file_name_(output_file_name),
    ppm_(ppm) {}

//sort feature 2d by intensities
FeaturePrsmPtr findMatchFeature(FeaturePrsmPtr feature, FeaturePrsmPtrVec &features, double ppm) {
  double cur_mass = feature->getMonoMass();
  double erro_tole = cur_mass * ppm;
  if (erro_tole < 0.01) {
    erro_tole = 0.01;
  }
  double cur_time = feature->getAlignRetentMiddle();
  for (size_t i = 0; i < features.size(); i++) {
    if (features[i] == nullptr) {
      continue;
    }
    double mass = features[i]->getMonoMass();
    if (abs(mass - cur_mass) > erro_tole) {
      continue;
    }

    double time = features[i]->getAlignRetentMiddle();
    if (abs(time-cur_time) < 0.02) {
      FeaturePrsmPtr result = features[i];
      features[i] = nullptr;
      return result;
    }
  }
  return nullptr;
}

bool featureUsed(FeaturePrsmPtr feature, FeaturePrsmPtrVec2D &table) {
  int cur_sample = feature->getSampleId();
  int feature_id = feature->getId();
  for (size_t i = 0; i < table.size(); i++) {
    FeaturePrsmPtr cur_feature= table[i][cur_sample];
    if (cur_feature != nullptr) {
      int cur_id = cur_feature->getId();
      if (feature_id == cur_id) {
        return true;
      }
    }
  }
  return false;
}

void getFeatureTable(FeaturePrsmPtrVec& all_features, FeaturePrsmPtrVec2D& features_2d, 
                     FeaturePrsmPtrVec2D& table, FeaturePrsmPtrVec &examples, double ppm) {
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
        cur_row[j] = findMatchFeature(all_features[i], features_2d[j],ppm);
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
  file << "Protein accession" << ","
      << "Protein description" << ","
      << "First residue" << ","
      << "Last residue" << ","
      << "Proteoform" << ","
      << "Precursor mass" << ",";

  for (int i = 0; i < sample_num; i++) {
    file << input_file_names_[i] << " Abundance" << ","
        << input_file_names_[i] << " Scan" << ","
        << input_file_names_[i] << " Retention time begin" << ","
        << input_file_names_[i] << " Retention time end" << ",";
  }
  file << std::endl;
  int cluster_num = table.size();
  for (int i = 0; i < cluster_num; i++) {
    FeaturePrsmPtr feature_ptr = examples[i];
    file << feature_ptr->getProtName() << ","
        << "\"" << feature_ptr->getProtDesc() << "\"" << ","
        << (feature_ptr->getFirstResidue() + 1) << ","
        << (feature_ptr->getLastResidue() + 1) << ","
        << feature_ptr->getProteoform() << ","
        << feature_ptr->getPrecMass() << ",";
    for (int j = 0; j < sample_num; j++) {
      FeaturePrsmPtr sample_feature = table[i][j];
      if (sample_feature == nullptr) {
        file << "," << "," << "," << ",";
      }
      else {
        file <<  std::setprecision(3) << std::scientific << sample_feature->getIntensity() << ","
            << std::fixed << sample_feature->getMs2Scan() << ","
            << sample_feature->getRetentBegin() << ","
            << sample_feature->getRetentEnd() << ",";
      }
    }
    file << std::endl;
  }
  file.close();
}

class Cell {
 public:
  Cell() {
    dist_ = -1;
    pre_ = -1;
  }
  double dist_ = -1;
  // 0: left, 1: up, 2, diag
  double pre_ = -1;
};

typedef std::shared_ptr<Cell> CellPtr;
typedef std::vector<CellPtr> CellPtrVec;
typedef std::vector<CellPtrVec> CellPtrVec2D;

void normalizeTime(FeaturePrsmPtrVec &features) {
  double max_time = 1;
  for (size_t i = 0; i < features.size(); i++) {
    if (features[i]->getRetentEnd() > max_time) {
      max_time = features[i]->getRetentEnd();
    }
  } 
  for (size_t i = 0; i < features.size(); i++) {
    double time_begin = features[i]->getRetentBegin()/max_time;
    features[i]->setAlignRetentBegin(time_begin);
    double time_end = features[i]->getRetentEnd()/max_time;
    features[i]->setAlignRetentEnd(time_end);
  }
}

inline FeaturePrsmPtrVec getTopFeatures(FeaturePrsmPtrVec &features) {
  std::sort(features.begin(), features.end(), Feature::cmpInteDec);
  FeaturePrsmPtrVec::const_iterator first = features.begin();
  FeaturePrsmPtrVec::const_iterator last = features.begin() + 1000;
  FeaturePrsmPtrVec result(first, last);
  return result;
}

inline CellPtrVec2D initCells(size_t a_size, size_t b_size, FeaturePrsmPtrVec &first) {
  CellPtrVec2D cell_2d;
  for (size_t i = 0; i < a_size + 1; i++) {
    CellPtrVec cells(b_size + 1);
    for (size_t j = 0; j < b_size + 1; j++) {
      cells[j] = std::make_shared<Cell>();
    }
    cell_2d.push_back(cells);
  }
  cell_2d[0][0]->dist_ = 0;
  for (size_t i = 1; i < a_size + 1; i++) {
    double inte = first[i-1]->getIntensity();
    cell_2d[i][0]->dist_ = cell_2d[i-1][0]->dist_ + std::log10(inte);
    cell_2d[i][0]->pre_ = 1; // up
  }
  for (size_t j = 1; j < b_size; j++) {
    cell_2d[0][j]->dist_ = 0;
    cell_2d[0][j]->pre_ = 0; // left
  }
  return cell_2d;
}

bool isMatch(FeaturePrsmPtr a, FeaturePrsmPtrVec &second, size_t b_index, double ppm) {
  double a_mass = a->getMonoMass();
  double erro_tole = a_mass * ppm;
  if (erro_tole < 0.01) {
    erro_tole = 0.01;
  }
  double b_time = second[b_index]->getAlignRetentMiddle();
  for (size_t i = b_index; i >= 1; i--) {
    double mass = second[i]->getMonoMass();
    double time = second[i]->getAlignRetentMiddle();
    if (b_time - time > 0.02) {
      break;
    }
    double error = std::abs(a_mass - mass);
    if (error <= erro_tole) {
      return true;
    }
  }
  for (size_t i = b_index; i < second.size(); i++) {
    double mass = second[i]->getMonoMass();
    double time = second[i]->getAlignRetentMiddle();
    if (time - b_time > 0.02) {
      break;
    }
    double error = std::abs(a_mass - mass);
    if (error <= erro_tole) {
      return true;
    }
  }
  return false;
}

void sampleAlignTime(FeaturePrsmPtrVec &first, FeaturePrsmPtrVec &second, 
                     double ppm, std::vector<std::pair<double,double>> &time_pairs) {
  FeaturePrsmPtrVec a = getTopFeatures(first);
  std::sort(a.begin(), a.end(), Feature::cmpRetentInc);
  FeaturePrsmPtrVec b = getTopFeatures(second);
  std::sort(b.begin(), b.end(), Feature::cmpRetentInc);
  size_t a_size = a.size();
  size_t b_size = b.size();
  LOG_DEBUG("a size " << a_size << " b size " << b_size);
  CellPtrVec2D cell_2d = initCells(a_size, b_size, first);
  LOG_DEBUG("cell inited ");
  //alignment
  for (size_t i = 1; i < a_size + 1; i++) {
    LOG_DEBUG("alignment  " << i);
    double inte = first[i-1]->getIntensity();
    double log_inte =  std::log10(inte);
    for (size_t j = 1; j < b_size + 1; j++) {
      double left_dist = cell_2d[i][j-1]->dist_;
      double up_dist = cell_2d[i-1][j]->dist_ + log_inte;
      double diag_dist = cell_2d[i-1][j-1]->dist_;
      if (!isMatch(first[i-1], second, j-1, ppm)) {
        diag_dist += log_inte;
      }
      if (diag_dist <= left_dist && diag_dist <= up_dist) {
        cell_2d[i][j]->dist_ = diag_dist;
        cell_2d[i][j]->pre_ = 2; // diag
      }
      else if (up_dist <= left_dist) {
        cell_2d[i][j]->dist_ = up_dist;
        cell_2d[i][j]->pre_ = 1; // up
      }
      else {
        cell_2d[i][j]->dist_ = left_dist;
        cell_2d[i][j]->pre_ = 0; // left
      }
    }
  }
  LOG_DEBUG("alignment finished ");
  //backtracking
  size_t i = a_size;
  size_t j = b_size;
  std::pair<double,double> last_pair(1,1);
  time_pairs.push_back(last_pair);
  while (i > 0 && j > 0) {
    double a_time = first[i-1]->getAlignRetentMiddle();
    double b_time = second[i-1]->getAlignRetentMiddle();
    std::pair<double, double>time_pair(a_time,b_time);
    time_pairs.push_back(time_pair);
    int pre = cell_2d[i][j]->pre_;
    if (pre == 2) {// diag
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
  LOG_DEBUG("backtrack finished ");
  std::reverse(time_pairs.begin(), time_pairs.end());
}

double findAlignTime(double time, std::vector<std::pair<double,double>> &time_pairs) {
  for (size_t i = 0; i < time_pairs.size() - 1; i++) {
    double a = time_pairs[i].second;
    double b = time_pairs[i+1].second;
    if (time >= a && time <= b) {
      double align_a = time_pairs[i].first;
      double align_b = time_pairs[i+1].first;
      if (a == b) {
        return align_a;
      }
      else {
        return align_a + (align_b - align_a) * (time-a)/(b-a);
      }
    }
  }
  LOG_ERROR("Can not find align time!");
  return -1;
}

void setAlignTime(FeaturePrsmPtrVec &features, 
                  std::vector<std::pair<double,double>> &time_pairs) {
  for (size_t i = 0; i < features.size(); i++) {
    double time_begin = findAlignTime(features[i]->getAlignRetentBegin(), time_pairs);
    features[i]->setAlignRetentBegin(time_begin);
    double time_end = findAlignTime(features[i]->getAlignRetentEnd(), time_pairs);
    features[i]->setAlignRetentEnd(time_end);
  }
}

void FeatureSampleMerge::process() {
  LOG_DEBUG("start process ");
  FeaturePrsmPtrVec2D features_2d; 
  FeaturePrsmPtrVec all_features; 
  size_t sample_num = input_file_names_.size();
  for (size_t k = 0; k < sample_num; k++) {
    std::string input_file_name = input_file_names_[k];
    FeaturePrsmReader reader(input_file_name);
    FeaturePrsmPtrVec features = reader.readAllFeatures();
    reader.close();
    std::sort(features.begin(), features.end(), Feature::cmpRetentInc);
    for (size_t i = 0; i < features.size(); i++) {
      features[i]->setSampleId(k);
    }
    normalizeTime(features);
    if (k > 0) {
      std::vector<std::pair<double,double>> time_pairs;
      sampleAlignTime(features_2d[0], features, ppm_, time_pairs);
      //setAlignTime(features, time_pairs);
    }
    features_2d.push_back(features);
    all_features.insert(all_features.end(), features.begin(), features.end());
  }
  std::sort(all_features.begin(), all_features.end(), Feature::cmpInteDec);
  for (size_t k = 0; k < sample_num; k++) {
    std::sort(features_2d[k].begin(), features_2d[k].end(), Feature::cmpInteDec);
  }
  FeaturePrsmPtrVec2D table; 
  FeaturePrsmPtrVec examples; 
  getFeatureTable(all_features, features_2d, table, examples, ppm_);
  outputTable(table, examples, sample_num);
}

}  // namespace prot


