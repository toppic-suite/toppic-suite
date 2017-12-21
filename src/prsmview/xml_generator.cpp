//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include "base/file_util.hpp"
#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_cluster.hpp"
#include "prsmview/anno_prsm.hpp"
#include "prsmview/anno_view.hpp"
#include "prsmview/xml_generator.hpp"
#include "spec/extend_ms_factory.hpp"

namespace prot {

XmlGenerator::XmlGenerator(PrsmParaPtr prsm_para_ptr,
                           const std::string &exec_dir,
                           const std::string &input_file_ext,
                           const std::string &fname_suffix) {
  input_file_ext_ = input_file_ext;
  mng_ptr_ = std::make_shared<PrsmViewMng>(prsm_para_ptr, exec_dir, fname_suffix);
  anno_view_ptr_ = std::make_shared<AnnoView>();
  fasta_reader_ptr_ = std::make_shared<FastaIndexReader>(prsm_para_ptr->getSearchDbFileName());
  writer_block_size_ = 300;
}

void XmlGenerator::outputPrsms() {
  std::set<int> species_id_set;
  std::set<int> prot_id_set;
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;

  PrsmReader prsm_reader(input_file_name);

  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                             mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getActivationPtr(),
                          mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getSkipList());
  SpectrumSetPtr spec_set_ptr;
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();

  size_t cnt = 0;
  size_t idx = 0;
  while ((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0])!= nullptr) {
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        species_id_set.insert(prsm_ptr->getProteoformPtr()->getProteoClusterId());
        prot_id_set.insert(prsm_ptr->getProteoformPtr()->getProtId());

        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        deconv_ms_vec2d_.push_back(deconv_ms_ptr_vec);
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        extend_ms_vec2d_.push_back(extend_ms_ptr_vec);
        spec_id_extend_ms_map_[spec_id] = idx;
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() +
            "prsms" + file_util::getFileSeparator() + "prsm" +
            string_util::convertToString(prsm_ptr->getPrsmId()) + ".xml";
        XmlWriter writer(file_name, "");
        writer.write(geneAnnoPrsm(writer.getDoc(), prsm_ptr, mng_ptr_));
        writer.close();

        std::vector<std::string> file_info;
        file_info.push_back(file_name);
        file_info.push_back(mng_ptr_->executive_dir_ + file_util::getFileSeparator() + "toppic_resources"
                            + file_util::getFileSeparator() + "xsl" + file_util::getFileSeparator() + "prsm.xsl");
        file_info.push_back(mng_ptr_->html_path_+ file_util::getFileSeparator() + "prsms" + file_util::getFileSeparator()
                            + "prsm" + string_util::convertToString(prsm_ptr->getPrsmId()) + ".html");
        anno_view_ptr_->file_list_.push_back(file_info); 

        cnt++;
        idx++;
        std::cout << std::flush << "Generating xml files - processing " << cnt << " PrSMs.\r";
        prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }
    }
  }

  std::cout << std::endl;

  prsm_reader.close();
  sp_reader.close();

  std::copy(species_id_set.begin(), species_id_set.end(), std::back_inserter(species_ids_));
  std::sort(species_ids_.begin(), species_ids_.end());
  std::copy(prot_id_set.begin(), prot_id_set.end(), std::back_inserter(prot_ids_));
  std::sort(prot_ids_.begin(), prot_ids_.end());
}

void XmlGenerator::outputProteoforms(){
  LOG_DEBUG("species id size " << species_ids_.size());
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();

  size_t cnt = 0;
  for (size_t i = 0; i < species_ids_.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating xml files - processing " << cnt << " proteoforms.\r";
    std::string input_file_name = file_util::basename(spectrum_file_name) + ".SPECIES_" + std::to_string(species_ids_[i]);
    PrsmPtrVec select_prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    for (size_t k = 0; k < select_prsm_ptrs.size(); k++) {
      select_prsm_ptrs[k]->setDeconvMsPtrVec(deconv_ms_vec2d_[spec_id_extend_ms_map_[select_prsm_ptrs[k]->getSpectrumId()]]);
      select_prsm_ptrs[k]->setRefineMsVec(extend_ms_vec2d_[spec_id_extend_ms_map_[select_prsm_ptrs[k]->getSpectrumId()]]);
    }
    if (select_prsm_ptrs.size() > 0) {
      std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() + "proteoforms" 
          + file_util::getFileSeparator() + "proteoform" + string_util::convertToString(species_ids_[i]) + ".xml";
      XmlWriter writer(file_name, "");
      std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(), Prsm::cmpEValueInc);
      writer.write(proteoformToXml(writer.getDoc(), select_prsm_ptrs, mng_ptr_));
      writer.close();
      LOG_DEBUG("output proteoform completed " << i);

      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_ptr_->executive_dir_ + file_util::getFileSeparator() + 
                          "toppic_resources" + file_util::getFileSeparator() + "xsl" + file_util::getFileSeparator() + "proteoform.xsl");
      file_info.push_back(mng_ptr_->html_path_+ file_util::getFileSeparator()+ "proteoforms" + file_util::getFileSeparator() 
                          + "proteoform"+string_util::convertToString(species_ids_[i])+".html");
      anno_view_ptr_->file_list_.push_back(file_info);
    }
  }
  std::cout << std::endl;
}

void XmlGenerator::outputProteins() {
  LOG_DEBUG("protein id size " << prot_ids_.size());
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  size_t cnt = 0;
  for (size_t i = 0; i < prot_ids_.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating xml files - processing " << std::ceil(cnt / 3.0) << " proteins.\r";
    int prot_id = prot_ids_[i]; 
    std::string input_file_name = file_util::basename(spectrum_file_name) + ".PROT_" + std::to_string(prot_id);
    PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    for (size_t k = 0; k < prsm_ptrs.size(); k++) {
      prsm_ptrs[k]->setDeconvMsPtrVec(deconv_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
      prsm_ptrs[k]->setRefineMsVec(extend_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
    }
    std::vector<int> species = prsm_util::getSpeciesIds(prsm_ptrs);
    if (!species.empty()) {
      std::string file_name = mng_ptr_->xml_path_ + file_util::getFileSeparator() + "proteins" 
          + file_util::getFileSeparator() + "protein" + string_util::convertToString(prot_id) + ".xml";
      XmlWriter writer(file_name,"");
      writer.write(proteinToXml(writer.getDoc(), prsm_ptrs, prot_id, species, mng_ptr_));
      writer.close();
      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_ptr_->executive_dir_ + file_util::getFileSeparator() + "toppic_resources" 
                          + file_util::getFileSeparator() + "xsl" + file_util::getFileSeparator() + "protein.xsl");
      file_info.push_back(mng_ptr_->html_path_+ file_util::getFileSeparator() + "proteins" + file_util::getFileSeparator() 
                          + "protein"+string_util::convertToString(prot_id)+".html");
      anno_view_ptr_->file_list_.push_back(file_info);
    }
  }
}

void XmlGenerator::outputAllProteins() {
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();

  PrsmPtrVec best_prsm_vec(prot_ids_.size());
  size_t cnt = prot_ids_.size();
  for (size_t i = 0; i < prot_ids_.size(); i++) {
    cnt++;
    std::cout << std::flush << "Generating xml files - processing " << std::ceil(cnt / 3.0) << " proteins.\r";
    std::string input_file_name = file_util::basename(spectrum_file_name) + ".PROT_" + std::to_string(prot_ids_[i]);
    PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    for (size_t k = 0; k < prsm_ptrs.size(); k++) {
      prsm_ptrs[k]->setDeconvMsPtrVec(deconv_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
      prsm_ptrs[k]->setRefineMsVec(extend_ms_vec2d_[spec_id_extend_ms_map_[prsm_ptrs[k]->getSpectrumId()]]);
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
    std::cout << std::flush << "Generating xml files - processing " << std::ceil(cnt / 3.0) << " proteins.\r";
    int prot_id = best_prsm_vec[i]->getProteoformPtr()->getProtId();
    std::string input_file_name = file_util::basename(spectrum_file_name) + ".PROT_" + std::to_string(prot_id);
    PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);
    std::vector<int> species = prsm_util::getSpeciesIds(prsm_ptrs);
    if (!species.empty()) {
      prot_elements->appendChild(proteinToXml(writer.getDoc(), prsm_ptrs, prot_id, species, mng_ptr_, false));
    }
  }
  std::cout << std::endl;
  writer.write(prot_elements);
  writer.close();
  std::vector<std::string> file_info;
  file_info.push_back(file_name);
  file_info.push_back(mng_ptr_->executive_dir_ + file_util::getFileSeparator() +
                      "toppic_resources" + file_util::getFileSeparator() + "xsl" + file_util::getFileSeparator() +"proteins.xsl");
  file_info.push_back(mng_ptr_->html_path_+ file_util::getFileSeparator() + "proteins.html");
  anno_view_ptr_->file_list_.push_back(file_info);
}

void XmlGenerator::outputFileList() {
  std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() + "files.xml";
  XmlWriter writer(file_name, "");
  writer.write(anno_view_ptr_->geneFileList(writer.getDoc()));
  writer.close();
}

void XmlGenerator::splitBySpeciesId() {
  std::cout << "Generating xml files - preprocessing " 
      << species_ids_.size() << " proteoforms." << std::endl;
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;

  std::vector<int> sub_species_ids_;
  std::vector<int>::iterator from = species_ids_.begin();
  bool flag = true;

  do {
    if (species_ids_.end() - from > writer_block_size_) {
      sub_species_ids_.assign(from, from + writer_block_size_);
      from += writer_block_size_;
    } else {
      sub_species_ids_.assign(from, species_ids_.end());
      flag = false;
    }
    std::map<int, size_t> species_id_map;
    std::vector<PrsmXmlWriterPtr> prsm_writer_vec;
    for (size_t i = 0; i < sub_species_ids_.size(); i++) {
      species_id_map[sub_species_ids_[i]] = i;
      prsm_writer_vec.push_back(std::make_shared<prot::PrsmXmlWriter>(file_util::basename(spectrum_file_name) + ".SPECIES_" + std::to_string(sub_species_ids_[i])));
    }

    PrsmReader prsm_reader(input_file_name);

    PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                               mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

    while (prsm_ptr != nullptr) {
      if (species_id_map.find(prsm_ptr->getProteoformPtr()->getProteoClusterId()) != species_id_map.end()) {
        prsm_writer_vec[species_id_map[prsm_ptr->getProteoformPtr()->getProteoClusterId()]]->write(prsm_ptr);
      }
      prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                         mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    }

    prsm_reader.close();

    for (size_t i = 0; i < prsm_writer_vec.size(); i++) {
      prsm_writer_vec[i]->close();
    }
  } while (flag);
}

void XmlGenerator::splitByProtId() {
  std::cout << "Generating xml files - preprocessing " 
      << prot_ids_.size() << " proteins." << std::endl;
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;

  std::vector<int> sub_prot_ids_;
  std::vector<int>::iterator from = prot_ids_.begin();

  bool flag = true;
  do {
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
      prsm_writer_vec.push_back(std::make_shared<prot::PrsmXmlWriter>(file_util::basename(spectrum_file_name) + ".PROT_" + std::to_string(sub_prot_ids_[i])));
    }

    PrsmReader prsm_reader(input_file_name);
    PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                               mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

    while (prsm_ptr != nullptr) {
      if (prot_id_map.find(prsm_ptr->getProteoformPtr()->getProtId()) != prot_id_map.end()) {
        prsm_writer_vec[prot_id_map[prsm_ptr->getProteoformPtr()->getProtId()]]->write(prsm_ptr);
      }
      prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                         mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    }

    prsm_reader.close();

    for (size_t i = 0; i < prsm_writer_vec.size(); i++) {
      prsm_writer_vec[i]->close();
    }  
  } while (flag);
}

void XmlGenerator::process() {
  LOG_DEBUG("process start");
  file_util::createFolder(mng_ptr_->xml_path_ + file_util::getFileSeparator() + "proteoforms");
  file_util::createFolder(mng_ptr_->xml_path_ + file_util::getFileSeparator() + "prsms");
  file_util::createFolder(mng_ptr_->xml_path_ + file_util::getFileSeparator() + "proteins");
  LOG_DEBUG("fold created");

  outputPrsms();
  splitBySpeciesId();
  outputProteoforms();

  splitByProtId();
  outputProteins();
  outputAllProteins();

  outputFileList();
}

}  // namespace prot
