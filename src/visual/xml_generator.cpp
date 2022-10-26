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

#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_cluster.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "visual/anno_prsm.hpp"
#include "visual/anno_file_list.hpp"
#include "visual/anno_xml_util.hpp"
#include "visual/xml_generator.hpp"

namespace toppic {

XmlGenerator::XmlGenerator(PrsmParaPtr prsm_para_ptr,
                           const std::string &resource_dir,
                           const std::string &input_file_ext,
                           const std::string &fname_suffix) {
  input_file_ext_ = input_file_ext;
  mng_ptr_ = std::make_shared<PrsmViewMng>(prsm_para_ptr, resource_dir, fname_suffix);
  anno_file_list_ptr_ = std::make_shared<AnnoFileList>();
  fasta_reader_ptr_ = std::make_shared<FastaIndexReader>(prsm_para_ptr->getOriDbName() + "_idx" + file_util::getFileSeparator() + prsm_para_ptr->getSearchDbFileName());
  writer_block_size_ = 300;
}

void XmlGenerator::outputPrsms() {
  std::set<int> cluster_id_set;
  std::set<int> prot_id_set;
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;

  PrsmReader prsm_reader(input_file_name);

  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                             mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  SimpleMsAlignReaderPtr ms_reader_ptr = std::make_shared<SimpleMsAlignReader>(sp_file_name, 
                                                                               group_spec_num,
                                                                               sp_para_ptr->getActivationPtr());

  SpectrumSetPtr spec_set_ptr;
  size_t cnt = 0;
  size_t idx = 0;
  while ((spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(ms_reader_ptr, sp_para_ptr))!= nullptr) {
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        cluster_id_set.insert(prsm_ptr->getProteoformPtr()->getProteoClusterId());
        prot_id_set.insert(prsm_ptr->getProteoformPtr()->getProtId());

        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        deconv_ms_vec2d_.push_back(deconv_ms_ptr_vec);
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        extend_ms_vec2d_.push_back(extend_ms_ptr_vec);
        spec_id_extend_ms_map_[spec_id] = idx;
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() +
            "prsms" + file_util::getFileSeparator() + "prsm" +
            str_util::toString(prsm_ptr->getPrsmId()) + ".xml";
        XmlWriter writer(file_name, "");
        writer.write(anno_prsm::geneAnnoPrsm(writer.getDoc(), prsm_ptr, mng_ptr_));
        writer.close();

        std::vector<std::string> file_info;
        file_info.push_back(file_name);
        file_info.push_back(mng_ptr_->html_path_ 
                              + file_util::getFileSeparator() + "data_js" 
                              + file_util::getFileSeparator() + "prsms" 
                              + file_util::getFileSeparator()
                              + "prsm" + str_util::toString(prsm_ptr->getPrsmId()) + ".js");
        anno_file_list_ptr_->file_list_.push_back(file_info); 
    
        cnt++;
        idx++;
        std::cout << std::flush << "Generating XML files - processing " << cnt << " single PrSMs.\r";
        prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_, 
                                           mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }
    }
  }

  std::cout << std::endl;

  prsm_reader.close();
  std::copy(cluster_id_set.begin(), cluster_id_set.end(), std::back_inserter(cluster_ids_));
  std::sort(cluster_ids_.begin(), cluster_ids_.end());
  std::copy(prot_id_set.begin(), prot_id_set.end(), std::back_inserter(prot_ids_));
  std::sort(prot_ids_.begin(), prot_ids_.end());
}

void XmlGenerator::outputAllPrsms() {
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx" + file_util::getFileSeparator() + mng_ptr_->prsm_para_ptr_->getSearchDbFileName();

  ModPtrVec fix_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();

  size_t cnt = 0;

  PrsmPtrVec prsm_ptrs 
      = prsm_reader_util::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
  for (size_t k = 0; k < prsm_ptrs.size(); k++) {
    prsm_ptrs[k]->setDeconvMsPtrVec(
        deconv_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
    prsm_ptrs[k]->setRefineMsVec(
        extend_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
  }
  std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() + "prsms.xml";
  XmlWriter writer(file_name, "");

  xercesc::DOMElement* prot_elements = writer.getDoc()->createElement("prsms");

  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating XML files - processing " << cnt << " PrSMs for a combined file.\r";
    prot_elements->appendChild(anno_prsm::geneAnnoPrsmBrief(writer.getDoc(), prsm_ptrs[i], mng_ptr_, true, false));                                                      
  }
  std::cout << std::endl;
  writer.write(prot_elements);
  writer.close();

  std::vector<std::string> file_info;
  file_info.push_back(file_name);
  file_info.push_back(mng_ptr_->html_path_+ file_util::getFileSeparator() + "data_js"
                        + file_util::getFileSeparator() + "prsms.js");
  anno_file_list_ptr_->file_list_.push_back(file_info);
}
void XmlGenerator::splitByProteoformId() {
  std::cout << "Generating XML files - preprocessing " 
      << cluster_ids_.size() << " proteoforms." << std::endl;
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;

  std::vector<int> sub_cluster_ids_;
  std::vector<int>::iterator from = cluster_ids_.begin();
  bool flag = true;

  while (flag) {
    if (cluster_ids_.end() - from > writer_block_size_) {
      sub_cluster_ids_.assign(from, from + writer_block_size_);
      from += writer_block_size_;
    } else {
      sub_cluster_ids_.assign(from, cluster_ids_.end());
      flag = false;
    }
    std::map<int, size_t> cluster_id_map;
    std::vector<PrsmXmlWriterPtr> prsm_writer_vec;
    for (size_t i = 0; i < sub_cluster_ids_.size(); i++) {
      cluster_id_map[sub_cluster_ids_[i]] = i;
      std::string file_name = file_util::basename(spectrum_file_name) 
          + ".proteoform_" + str_util::toString(sub_cluster_ids_[i]);
      PrsmXmlWriterPtr writer_ptr = std::make_shared<toppic::PrsmXmlWriter>(file_name);
      prsm_writer_vec.push_back(writer_ptr);
    }

    PrsmReader prsm_reader(input_file_name);

    PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                               mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

    while (prsm_ptr != nullptr) {
      int prsm_cluster_id = prsm_ptr->getProteoformPtr()->getProteoClusterId();
      if (cluster_id_map.find(prsm_cluster_id) != cluster_id_map.end()) {
        prsm_writer_vec[cluster_id_map[prsm_cluster_id]]->write(prsm_ptr);
      }
      prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                         mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    }

    prsm_reader.close();

    for (size_t i = 0; i < prsm_writer_vec.size(); i++) {
      prsm_writer_vec[i]->close();
    }
  } 
}

void XmlGenerator::outputProteoforms(){
  LOG_DEBUG("prsm cluster id size " << cluster_ids_.size());
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  //std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx" + file_util::getFileSeparator() + mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();

  size_t cnt = 0;
  for (size_t i = 0; i < cluster_ids_.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating XML files - processing " << cnt << " proteoforms.\r";
    std::string input_file_name = file_util::basename(spectrum_file_name) 
        + ".proteoform_" + str_util::toString(cluster_ids_[i]);
    PrsmPtrVec select_prsm_ptrs 
        = prsm_reader_util::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    for (size_t k = 0; k < select_prsm_ptrs.size(); k++) {
      select_prsm_ptrs[k]->setDeconvMsPtrVec(
          deconv_ms_vec2d_[spec_id_extend_ms_map_[select_prsm_ptrs[k]->getSpectrumId()]]);
      select_prsm_ptrs[k]->setRefineMsVec(
          extend_ms_vec2d_[spec_id_extend_ms_map_[select_prsm_ptrs[k]->getSpectrumId()]]);
    }
    if (select_prsm_ptrs.size() > 0) {
      std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() 
          + "proteoforms" + file_util::getFileSeparator() 
          + "proteoform" + str_util::toString(cluster_ids_[i]) + ".xml";
      XmlWriter writer(file_name, "");
      std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(), Prsm::cmpEValueInc);
      bool detail = true; 
      bool add_ms = true;
      writer.write(anno_xml_util::geneXmlForProteoform(writer.getDoc(), select_prsm_ptrs, 
                                                       mng_ptr_, detail, add_ms));
      writer.close();
      LOG_DEBUG("output proteoform completed " << i);
        std::vector<std::string> file_info;
        file_info.push_back(file_name);
        file_info.push_back(mng_ptr_->html_path_ 
                            + file_util::getFileSeparator() + "data_js" 
                            + file_util::getFileSeparator() + "proteoforms" 
                            + file_util::getFileSeparator()
                            + "proteoform"+str_util::toString(cluster_ids_[i])+".js");
        anno_file_list_ptr_->file_list_.push_back(file_info);
    }
  }
  std::cout << std::endl;
}

void XmlGenerator::splitByProtId() {
  std::cout << "Generating XML files - preprocessing " 
      << prot_ids_.size() << " proteins." << std::endl;
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;

  std::vector<int> sub_prot_ids_;
  std::vector<int>::iterator from = prot_ids_.begin();

  bool flag = true;
  while (flag) {
    if (prot_ids_.end() - from > writer_block_size_) {
      sub_prot_ids_.assign(from, from + writer_block_size_);
      from += writer_block_size_;
    } else {
      sub_prot_ids_.assign(from, prot_ids_.end());
      flag = false;
    }
    std::map<int, size_t> prot_id_map;
    std::vector<PrsmXmlWriterPtr> prsm_writer_vec;
    for (size_t i = 0; i < sub_prot_ids_.size(); i++) {
      prot_id_map[sub_prot_ids_[i]] = i;
      std::string file_name = file_util::basename(spectrum_file_name) 
          + ".prot_" + str_util::toString(sub_prot_ids_[i]);
      PrsmXmlWriterPtr writer_ptr = std::make_shared<toppic::PrsmXmlWriter>(file_name);
      prsm_writer_vec.push_back(writer_ptr);
    }

    PrsmReader prsm_reader(input_file_name);
    PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                               mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

    while (prsm_ptr != nullptr) {
      int prot_id = prsm_ptr->getProteoformPtr()->getProtId();
      if (prot_id_map.find(prot_id) != prot_id_map.end()) {
        prsm_writer_vec[prot_id_map[prot_id]]->write(prsm_ptr);
      }
      prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                         mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    }

    prsm_reader.close();

    for (size_t i = 0; i < prsm_writer_vec.size(); i++) {
      prsm_writer_vec[i]->close();
    }  
  }
}

void XmlGenerator::outputProteins() {
  LOG_DEBUG("protein id size " << prot_ids_.size());
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  //std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx" + file_util::getFileSeparator() + mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  size_t cnt = 0;
  for (size_t i = 0; i < prot_ids_.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating XML files - processing " 
        << std::ceil(cnt / 3.0) << " proteins.\r";
    int prot_id = prot_ids_[i]; 
    std::string input_file_name = file_util::basename(spectrum_file_name) 
        + ".prot_" + str_util::toString(prot_id);
    PrsmPtrVec prsm_ptrs 
        = prsm_reader_util::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    for (size_t k = 0; k < prsm_ptrs.size(); k++) {
      prsm_ptrs[k]->setDeconvMsPtrVec(
          deconv_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
      prsm_ptrs[k]->setRefineMsVec(
          extend_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
    }
    std::vector<int> cluster = prsm_util::getClusterIds(prsm_ptrs);
    if (!cluster.empty()) {
      std::string file_name = mng_ptr_->xml_path_ + file_util::getFileSeparator() + "proteins" 
          + file_util::getFileSeparator() + "protein" + str_util::toString(prot_id) + ".xml";
      XmlWriterPtr writer = std::make_shared<XmlWriter>(file_name, "");
      anno_xml_util::writeProteinToXml(writer, prsm_ptrs, prot_id, cluster, mng_ptr_, true, false);
      writer->close();

      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_ptr_->html_path_ 
                            + file_util::getFileSeparator() + "data_js" 
                            + file_util::getFileSeparator() + "proteins" 
                            + file_util::getFileSeparator() 
                            + "protein"+str_util::toString(prot_id)+".js");
      anno_file_list_ptr_->file_list_.push_back(file_info);
    }
  }
}

void XmlGenerator::outputAllProteins() {
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  //std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx" + file_util::getFileSeparator() + mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();

  PrsmPtrVec best_prsm_vec(prot_ids_.size());
  size_t cnt = 0; 
  for (size_t i = 0; i < prot_ids_.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating XML files - processing " 
        << std::ceil(cnt / 2.0) << " proteins.\r";
    std::string input_file_name = file_util::basename(spectrum_file_name) 
        + ".prot_" + str_util::toString(prot_ids_[i]);
    PrsmPtrVec prsm_ptrs 
        = prsm_reader_util::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    for (size_t k = 0; k < prsm_ptrs.size(); k++) {
      prsm_ptrs[k]->setDeconvMsPtrVec(
          deconv_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
      prsm_ptrs[k]->setRefineMsVec(
          extend_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
    }
    std::sort(prsm_ptrs.begin(), prsm_ptrs.end(), Prsm::cmpEValueInc);
    best_prsm_vec[i] = prsm_ptrs[0];
  }
  std::sort(best_prsm_vec.begin(), best_prsm_vec.end(), Prsm::cmpEValueInc);
  std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() + "proteins.xml";
  XmlWriter writer(file_name, "protein_list");

  xercesc::DOMElement* prot_elements = writer.getDoc()->createElement("proteins");

  for (size_t i = 0; i < best_prsm_vec.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating XML files - processing " 
        << std::ceil(cnt / 2.0) << " proteins.\r";
    int prot_id = best_prsm_vec[i]->getProteoformPtr()->getProtId();
    std::string input_file_name = file_util::basename(spectrum_file_name) 
        + ".prot_" + str_util::toString(prot_id);
    PrsmPtrVec prsm_ptrs 
        = prsm_reader_util::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    std::vector<int> cluster = prsm_util::getClusterIds(prsm_ptrs);
    if (!cluster.empty()) {
      prot_elements->appendChild(
          anno_xml_util::geneXmlForProteinList(writer.getDoc(), prsm_ptrs, prot_id, cluster, 
                                               mng_ptr_, false));                                     
    }
  }
  std::cout << std::endl;
  writer.write(prot_elements);
  writer.close();

  std::vector<std::string> file_info;
  file_info.push_back(file_name);
  file_info.push_back(mng_ptr_->html_path_+ file_util::getFileSeparator() + "data_js"
                        + file_util::getFileSeparator() + "proteins.js");
  anno_file_list_ptr_->file_list_.push_back(file_info);
  
}

void XmlGenerator::outputFileList() {
  std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() + "files.xml";
  XmlWriter writer(file_name, "");
  writer.write(anno_file_list_ptr_->geneFileList(writer.getDoc()));
  writer.close();
}

void XmlGenerator::removeTempFiles() {
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  file_util::cleanTempFiles(sp_file_name, "proteoform_");
  file_util::cleanTempFiles(sp_file_name, "prot_");
}

void XmlGenerator::process() {
  LOG_DEBUG("process start");
  if (file_util::exists(mng_ptr_->xml_path_)) {
    LOG_WARN("The XML directory " << mng_ptr_->xml_path_ << " exists!");
    file_util::delDir(mng_ptr_->xml_path_);
  }
  file_util::createFolder(mng_ptr_->xml_path_ + file_util::getFileSeparator() + "proteoforms");
  file_util::createFolder(mng_ptr_->xml_path_ + file_util::getFileSeparator() + "prsms");
  file_util::createFolder(mng_ptr_->xml_path_ + file_util::getFileSeparator() + "proteins");
  LOG_DEBUG("fold created");

  outputPrsms();
  outputAllPrsms();
  splitByProteoformId();
  outputProteoforms();

  splitByProtId();
  outputProteins();
  outputAllProteins();

  outputFileList();
  removeTempFiles();
}

}  // namespace toppic
