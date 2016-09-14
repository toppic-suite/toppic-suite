// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include <set>
#include <boost/algorithm/string.hpp>

#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/fasta_reader.hpp"
#include "base/proteoform_factory.hpp"
#include "spec/peak.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_util.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_prsm.hpp"
#include "prsmview/anno_view.hpp"

namespace prot{

xercesc::DOMElement* AnnoView::geneFileList(XmlDOMDocument* xml_doc){
  xercesc::DOMElement* element = xml_doc->createElement("file_list");
  for(size_t i=0;i<file_list_.size();i++){
    xercesc::DOMElement* file = xml_doc->createElement("file");
    xml_doc->addElement(file, "xml", file_list_[i][0].c_str());
    xml_doc->addElement(file, "xsl", file_list_[i][1].c_str());
    xml_doc->addElement(file, "html", file_list_[i][2].c_str());
    element->appendChild(file);
  }
  return element;
}

std::vector<std::vector<std::string>> readViewXmlFiles(const std::string &file_name){
  std::vector<std::vector<std::string>> file_list;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if(parser){
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name.c_str());
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      int file_num = XmlDomUtil::getChildCount(root, "file");
      for (int i = 0; i < file_num; i++) {
        xercesc::DOMElement* file_element = XmlDomUtil::getChildElement(root, "file", i);
        std::vector<std::string> file_info;
        file_info.push_back(XmlDomUtil::getChildValue(file_element,"xml",0));
        file_info.push_back(XmlDomUtil::getChildValue(file_element,"xsl",0));
        file_info.push_back(XmlDomUtil::getChildValue(file_element,"html",0));
        file_list.push_back(file_info);
      }
    }
    delete doc;
  }
  return file_list;
}

xercesc::DOMElement* proteoformToXml(XmlDOMDocument* xml_doc, const PrsmPtrVec &prsm_ptrs, 
                                     PrsmViewMngPtr mng_ptr, bool detail){
  xercesc::DOMElement* proteoform_element = xml_doc->createElement("compatible_proteoform");
  std::string str=StringUtil::convertToString(prsm_ptrs[0]->getProteoformPtr()->getProtId());
  xml_doc->addElement(proteoform_element, "sequence_id", str.c_str());
  str=prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(proteoform_element, "sequence_name", str.c_str());
  str=prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(proteoform_element, "sequence_description", str.c_str());
  str=StringUtil::convertToString(prsm_ptrs[0]->getProteoformPtr()->getSpeciesId());
  xml_doc->addElement(proteoform_element, "proteoform_id", str.c_str());
  int count = prsm_ptrs.size();
  str=StringUtil::convertToString(count);
  xml_doc->addElement(proteoform_element, "prsm_number", str.c_str());
  for(size_t i=0;i<prsm_ptrs.size();i++){
    proteoform_element->appendChild(geneAnnoPrsm(xml_doc,prsm_ptrs[i], mng_ptr, detail));
  }
  return proteoform_element;
}


xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  FastaSeqPtr seq_ptr,
                                  int prot_id,
                                  const std::vector<int> &species_ids,
                                  PrsmViewMngPtr mng_ptr,
                                  bool detail){
  xercesc::DOMElement* prot_element = xml_doc->createElement("protein");
  std::string str=StringUtil::convertToString(prot_id);
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str=seq_ptr->getName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str=seq_ptr->getDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  int count = species_ids.size();
  str=StringUtil::convertToString(count);
  xml_doc->addElement(prot_element, "compatible_proteoform_number", str.c_str());
  for(size_t i=0;i<species_ids.size();i++){
    PrsmPtrVec select_prsm_ptrs = PrsmUtil::selectSpeciesPrsms(prsm_ptrs,species_ids[i]);
    std::sort(select_prsm_ptrs.begin(),select_prsm_ptrs.end(),Prsm::cmpEValueInc);
    prot_element->appendChild(proteoformToXml(xml_doc,select_prsm_ptrs, mng_ptr, detail));
  }
  return prot_element;
}

PrsmPtr getBestEValuePrsmPtr (std::string &seq_name, const PrsmPtrVec &prsm_ptrs) {
  PrsmPtr best_ptr(nullptr);
  double best_evalue = std::numeric_limits<double>::max();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getSeqName() == seq_name && 
        prsm_ptrs[i]->getEValue() < best_evalue) {
      best_evalue = prsm_ptrs[i]->getEValue();
      best_ptr = prsm_ptrs[i];
    }
  }
  return best_ptr;
}

inline bool evalueCompare(const std::pair<FastaSeqPtr, double> &a, const std::pair<FastaSeqPtr, double> &b) {
  return a.second < b.second;
}


xercesc::DOMElement* allProteinToXml(XmlDOMDocument* xml_doc,
                                     const PrsmPtrVec &prsm_ptrs,
                                     PrsmViewMngPtr mng_ptr){
  xercesc::DOMElement* prot_elements = xml_doc->createElement("proteins");
  // sort 
  FastaReader reader(mng_ptr->prsm_para_ptr_->getSearchDbFileName());
  FastaSeqPtr seq_ptr = reader.getNextSeq();
  std::vector<std::pair<FastaSeqPtr, double>> seq_evalues;

  while (seq_ptr != nullptr) {
    mng_ptr->cnt_++;
    std::cout << std::flush << "Generating xml files - processing " << mng_ptr->cnt_ << " of " << mng_ptr->num_files_ << " files.\r";
    std::string seq_name = seq_ptr->getName();
    PrsmPtr best_ptr = getBestEValuePrsmPtr (seq_name, prsm_ptrs);
    if (best_ptr != nullptr) {
      std::pair<FastaSeqPtr, double> cur_seq_evalue(seq_ptr, best_ptr->getEValue());
      seq_evalues.push_back(cur_seq_evalue);
    }
    seq_ptr = reader.getNextSeq();
  }
  std::sort(seq_evalues.begin(), seq_evalues.end(), evalueCompare);

  for(size_t i=0;i<seq_evalues.size();i++){
    std::string seq_name = seq_evalues[i].first->getName();
    std::vector<int> species_ids = PrsmUtil::getSpeciesIds(prsm_ptrs,seq_name);
    int prot_id = PrsmUtil::getProteinId(prsm_ptrs, seq_name);
    if(species_ids.size()>0){
      prot_elements->appendChild(proteinToXml(xml_doc,prsm_ptrs,seq_evalues[i].first, 
                                              prot_id, species_ids, mng_ptr, false));
    }
  }
  return prot_elements;
}

}

