//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include <set>
#include <algorithm>
#include <map>

#include "common/util/logger.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/feature/spec_feature_reader.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_util.hpp"

namespace toppic {

namespace prsm_util {

void setValueStr(std::vector<std::string> &str_vec, const std::string &property, std::string val) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      std::string line = str_vec[i];
      int start = line.find(">");
      int end = line.find("<", start);
      std::string start_tag = line.substr(0, start + 1);
      std::string end_tag = line.substr(end);
      str_vec[i] = start_tag + val + end_tag;
    }
  }
};

std::string getValueStr(std::string line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  return num_str;
}

std::string getXmlLine(const std::vector<std::string> &str_vec,
                       const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return str_vec[i];
    }
  }
  return "";
}

std::vector<std::string> getXmlLineVec(const std::vector<std::string> &str_vec,
                                       const std::string &property) {
  std::vector<std::string> res;
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      res.push_back(str_vec[i]);
    }
  }
  return res;
}


PrsmPtrVec selectClusterPrsms(const PrsmPtrVec &prsm_ptrs, int cluster_id) {
  PrsmPtrVec select_prsm_ptrs;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (cluster_id == prsm_ptrs[i]->getProteoformPtr()->getProteoClusterId()) {
      select_prsm_ptrs.push_back(prsm_ptrs[i]);
    }
  }
  return select_prsm_ptrs;
}

std::vector<int> getClusterIds(const PrsmPtrVec &prsm_ptrs, std::string &seq_name) {
  std::set<int> cluster_id_set;
  std::vector<int> cluster_ids;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getSeqName() == seq_name)
      cluster_id_set.insert(prsm_ptrs[i]->getProteoformPtr()->getProteoClusterId());
  }
  std::copy(cluster_id_set.begin(), cluster_id_set.end(), std::back_inserter(cluster_ids));
  std::sort(cluster_ids.begin(), cluster_ids.end());
  return cluster_ids;
}

int getProteinId(const PrsmPtrVec &prsm_ptrs, const std::string &seq_name) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getSeqName() == seq_name) {
      return prsm_ptrs[i]->getProteoformPtr()->getProtId();
    }
  }
  return -1;
}

std::vector<int> getClusterIds(const PrsmPtrVec &prsm_ptrs) {
  std::set<int> cluster_id_set;
  std::vector<int> cluster_ids;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    cluster_id_set.insert(prsm_ptrs[i]->getProteoformPtr()->getProteoClusterId());
  }
  std::copy(cluster_id_set.begin(), cluster_id_set.end(), std::back_inserter(cluster_ids));
  std::sort(cluster_ids.begin(), cluster_ids.end());
  return cluster_ids;
}

bool isMatchMs(PrsmPtr prsm_ptr, MsHeaderPtr header_ptr) {
  int id = header_ptr->getSpecId();
  std::string scan = header_ptr->getScansString();
  int prec_id = header_ptr->getFirstPrecId();
  if (prsm_ptr->getSpectrumId() == id && prsm_ptr->getPrecursorId() == prec_id) {
    if (prsm_ptr->getSpectrumScan() != scan) {
      LOG_ERROR("Error in Prsm! Spectrum id:" << prsm_ptr->getSpectrumId());
    }
    return true;
  } else {
    return false;
  }
}

void addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr) {
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  MsAlignReaderPtr ms_reader_ptr = std::make_shared<MsAlignReader>(prsm_para_ptr->getSpectrumFileName(),
                                                                   prsm_para_ptr->getGroupSpecNum(),
                                                                   sp_para_ptr->getActivationPtr());
  SpectrumSetPtr spec_set_ptr 
      = spectrum_set_factory::readNextSpectrumSetPtr(ms_reader_ptr, sp_para_ptr);
  // use prsm order information (ordered by spectrum id then prec id)
  int start_prsm = 0;
  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
      int spectrum_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getSpecId();
      int prec_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecId();
      LOG_DEBUG("spectrum id " << spectrum_id);
      for (size_t i = start_prsm; i < prsm_ptrs.size(); i++) {
        if (isMatchMs(prsm_ptrs[i], deconv_ms_ptr_vec[0]->getMsHeaderPtr())) {
          prsm_ptrs[i]->setDeconvMsPtrVec(deconv_ms_ptr_vec);
          double new_prec_mass = prsm_ptrs[i]->getAdjustedPrecMass();
          prsm_ptrs[i]->setRefineMsVec(
              extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec,
                                                   prsm_para_ptr->getSpParaPtr(),
                                                   new_prec_mass));
        }
        if ((spectrum_id == prsm_ptrs[i]->getSpectrumId()
             && prec_id < prsm_ptrs[i]->getPrecursorId()) ||
            spectrum_id < prsm_ptrs[i]->getSpectrumId()) {
          start_prsm = i;
          break;
        }
      }
    }
    spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(ms_reader_ptr, sp_para_ptr);
  }
}

void addFeatureInfoToPrsms(PrsmStrPtrVec &prsm_ptrs, const std::string & feature_file_name) {
  // read TopFD featuers
  SpecFeatureReader ft_reader(feature_file_name); 
  SpecFeaturePtrVec ms2_features = ft_reader.readAllFeatures();
  ft_reader.close();

  std::map<int,SpecFeaturePtr> feature_map;
  for (size_t i = 0; i < ms2_features.size(); i++) {
    feature_map[ms2_features[i]->getFracFeatureId()] =  ms2_features[i];
  }

  // make sure prsms sorted by spectrum id
  std::sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);

  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int spec_id = prsm_ptrs[i]->getSpectrumId();
    int frac_feature_id = prsm_ptrs[i]->getFracFeatureId();
    //LOG_ERROR(spec_id << " " << sample_feature_id);
    if (feature_map.find(frac_feature_id) != feature_map.end()) { 
      SpecFeaturePtr feature = feature_map.find(frac_feature_id)->second;
      if (feature != nullptr) {
        prsm_ptrs[i]->setFracFeatureInte(feature->getFracFeatureInte());
        prsm_ptrs[i]->setFracFeatureScore(feature->getFracFeatureScore());
        prsm_ptrs[i]->setFracFeatureApexTime(feature->getFracFeatureApexTime());
        prsm_ptrs[i]->setFracFeatureMinTime(feature->getFracFeatureMinTime());
        prsm_ptrs[i]->setFracFeatureMaxTime(feature->getFracFeatureMaxTime());
      }
      else {
        LOG_ERROR("Spectrum " << spec_id << " does not have a feature!");
      }
    }
  }
}

void removePrsmsWithoutFeature(PrsmStrPtrVec &prsm_ptrs, 
                               PrsmStrPtrVec &filtered_prsm_ptrs) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getFracFeatureId() >=0) {
      filtered_prsm_ptrs.push_back(prsm_ptrs[i]);
    }
  }
}


void mergePrsmFiles(const std::vector<std::string> & prsm_file_lst, 
                    int max_spec_num_per_file,
                    int max_feat_num_per_file,
                    const std::string & output_file) {
  PrsmXmlWriterPtr prsm_writer = std::make_shared<PrsmXmlWriter>(output_file);

  for (size_t i = 0; i < prsm_file_lst.size(); i++) {
    PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(prsm_file_lst[i]); 
    PrsmStrPtr prsm = prsm_reader->readOnePrsmStr();
    while (prsm != nullptr) {
      prsm->setSpectrumId(max_spec_num_per_file * i + prsm->getSpectrumId());
      prsm->setFracFeatureId(max_feat_num_per_file * i + prsm->getFracFeatureId());
      prsm_writer->write(prsm);
      prsm = prsm_reader->readOnePrsmStr();
    }
  }

  prsm_writer->close();
}

double compClusterInte(PrsmStrPtrVec prsm_list) {
  double inte = 0;
  std::set<int> feat_id_set;
  for (size_t i = 0; i < prsm_list.size(); i++) {
    PrsmStrPtr prsm_ptr = prsm_list[i];
    int feat_id = prsm_ptr->getFracFeatureId();
    if (feat_id_set.find(feat_id) == feat_id_set.end()) {
      feat_id_set.insert(feat_id);
      inte = inte + prsm_ptr->getFracFeatureInte();
    }
  }
  return inte;
}

}  // namespace prsm_util

}  // namespace toppic
