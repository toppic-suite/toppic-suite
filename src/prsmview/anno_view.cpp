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
                                     PrsmViewMngPtr mng_ptr){
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
    proteoform_element->appendChild(geneAnnoPrsm(xml_doc,prsm_ptrs[i], mng_ptr));
  }
  return proteoform_element;
}


xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  FastaSeqPtr seq_ptr,
                                  int prot_id,
                                  const std::vector<int> &species_ids,
                                  PrsmViewMngPtr mng_ptr){
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
    prot_element->appendChild(proteoformToXml(xml_doc,select_prsm_ptrs, mng_ptr));
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
      prot_elements->appendChild(proteinToXml(xml_doc,prsm_ptrs,seq_evalues[i].first, prot_id, species_ids, mng_ptr));
    }
  }
  return prot_elements;
}

}

/*
std::vector<xercesc::DOMElement*> modificationToXml(XmlDOMDocument* xml_doc,
        PrsmPtrVec & prsm_ptrs, const PrsmViewMngPtr & mng_ptr) {
    double ppo = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
    xercesc::DOMElement* mod_element;
    std::set<std::string> mod_set;

    for (size_t i = 0; i < prsm_ptrs.size(); i++) {
        double err = ppo * prsm_ptrs[i]->getOriPrecMass();
        ChangePtrVec change_vec = prsm_ptrs[i]->getProteoformPtr()->getUnexpectedChangePtrVec(err);
        for (size_t j = 0; j < change_vec.size(); j++) {
            if (change_vec[j]->getPtmPtr() != nullptr)
                mod_set.insert(change_vec[j]->getPtmPtr()->getName());
        }
    }

    std::vector<xercesc::DOMElement *> xml_vec;

    for (std::set<std::string>::iterator iter = mod_set.begin(); iter != mod_set.end(); iter++) {
        PrsmPtrVec select_prsm_ptrs;
        std::for_each(prsm_ptrs.begin(), prsm_ptrs.end(),
        [iter, &select_prsm_ptrs, ppo](PrsmPtr p) {
            double err = ppo * p->getOriPrecMass();
            ChangePtrVec change_vec = p->getProteoformPtr()->getUnexpectedChangePtrVec(err);
            for (size_t i = 0; i < change_vec.size(); i++) {
                if (change_vec[i]->getPtmPtr() != nullptr) {
                    if (change_vec[i]->getPtmPtr()->getName() == *iter) {
                        select_prsm_ptrs.push_back(p);
                        break;
                    }
                }
            }

        });

        mod_element = xml_doc->createElement("modification");
        xml_doc->addElement(mod_element, "modification_name", (*iter).c_str());
        std::sort(select_prsm_ptrs.begin(),select_prsm_ptrs.end(),prsmEValueUp);
        mod_element->appendChild(proteoformToXml(xml_doc,select_prsm_ptrs, mng_ptr));

        xml_vec.push_back(mod_element);
    }
    return xml_vec;
}
*/

/*
xercesc::DOMElement* allModificationToXml(XmlDOMDocument* xml_doc,
        PrsmPtrVec & prsm_ptrs, const PrsmViewMngPtr & mng_ptr) {

    double ppo = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
    xercesc::DOMElement* mod_elements = xml_doc->createElement("modifications");

    std::sort(prsm_ptrs.begin(), prsm_ptrs.end(),
    [](const PrsmPtr & lhs, const PrsmPtr & rhs) {
        return lhs->getProteoformPtr()->getDbResSeqPtr()->getId() < rhs->getProteoformPtr()->getDbResSeqPtr()->getId();
    });

    int prot_id = prsm_ptrs[0]->getProteoformPtr()->getDbResSeqPtr()->getId();

    PrsmPtrVec prsm_mod;
    std::vector<xercesc::DOMElement*> mods;
    for (size_t i = 0; i< prsm_ptrs.size(); i++) {

        double err = ppo * prsm_ptrs[i]->getOriPrecMass();
        if (prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getId() == prot_id) {
            if (prsm_ptrs[i]->getProteoformPtr()->getUnexpectedChangeNum(err) > 0)
                prsm_mod.push_back(prsm_ptrs[i]);
        } else {
            mods.clear();
            mods = modificationToXml(xml_doc, prsm_mod, mng_ptr);
            for (size_t i = 0; i < mods.size(); i++) {
                mod_elements->appendChild(mods[i]);
            }
            prsm_mod.clear();
            prsm_mod.push_back(prsm_ptrs[i]);
            prot_id = prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getId();
        }
    }

    mods = modificationToXml(xml_doc, prsm_mod, mng_ptr);

    for (size_t i = 0; i < mods.size(); i++) {
        mod_elements->appendChild(mods[i]);
    }
    return mod_elements;
}
*/

