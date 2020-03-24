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

#include <set>
#include <algorithm>
#include <map>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/feature/spec_feature_reader.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_util.hpp"

namespace toppic {

namespace prsm_util {

std::string getValueStr(std::string line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  // LOG_DEBUG(line << "  Num str: " << num_str);
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
  int id = header_ptr->getId();
  std::string scan = header_ptr->getScansString();
  int prec_id = header_ptr->getPrecId();
  if (prsm_ptr->getSpectrumId() == id && prsm_ptr->getPrecursorId() == prec_id) {
    if (prsm_ptr->getSpectrumScan() != scan) {
      LOG_ERROR("Error in Prsm.");
    }
    return true;
  } else {
    LOG_TRACE("prsm spectrum id " << prec_id << " ms spectrum id " << id);
    return false;
  }
}

void addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr) {
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  SimpleMsAlignReaderPtr ms_reader_ptr 
      = std::make_shared<SimpleMsAlignReader>(prsm_para_ptr->getSpectrumFileName(),
                                              prsm_para_ptr->getGroupSpecNum(),
                                              sp_para_ptr->getActivationPtr());
  SpectrumSetPtr spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(ms_reader_ptr, sp_para_ptr);
  // use prsm order information (ordered by spectrum id then prec id)
  int start_prsm = 0;
  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
      int spectrum_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getId();
      int prec_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecId();
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

void addFeatureIDToPrsms(PrsmStrPtrVec &prsm_ptrs, const std::string & feature_file_name) {
  // read TopFD featuers
  SpecFeatureReader ft_reader(feature_file_name); 
  SpecFeaturePtrVec ms2_features = ft_reader.readAllFeatures();
  ft_reader.close();

  std::map<int,SpecFeaturePtr> feature_map;
  for (size_t i = 0; i < ms2_features.size(); i++) {
    feature_map[ms2_features[i]->getSpecId()] =  ms2_features[i];
  }

  // make sure prsms sorted by spectrum id
  std::sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpSpectrumIdInc);

  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int spec_id = prsm_ptrs[i]->getSpectrumId();
    if (feature_map.find(spec_id) != feature_map.end()) { 
      SpecFeaturePtr feature = feature_map.find(spec_id)->second;
      if (feature != nullptr) {
        prsm_ptrs[i]->setPrecFeatureId(feature->getSampleFeatureId());
        prsm_ptrs[i]->setPrecFeatureInte(feature->getSampleFeatureInte());
        prsm_ptrs[i]->setFracFeatureScore(feature->getFracFeatureScore());
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
    if (prsm_ptrs[i]->getSampleFeatureId() >=0) {
      filtered_prsm_ptrs.push_back(prsm_ptrs[i]);
    }
  }
}


void mergePrsmFiles(const std::vector<std::string> & prsm_file_lst, int N, 
                    const std::string & output_file) {
  PrsmXmlWriterPtr prsm_writer = std::make_shared<PrsmXmlWriter>(output_file);

  for (size_t i = 0; i < prsm_file_lst.size(); i++) {
    PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(prsm_file_lst[i]); 
    PrsmStrPtr prsm = prsm_reader->readOnePrsmStr();
    while (prsm != nullptr) {
      prsm->setSpectrumId(N * i + prsm->getSpectrumId());

      prsm_writer->write(prsm);

      prsm = prsm_reader->readOnePrsmStr();
    }
  }

  prsm_writer->close();
}

}  // namespace prsm_util

}  // namespace toppic
