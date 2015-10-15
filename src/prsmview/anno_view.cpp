
#include <set>
#include <boost/algorithm/string.hpp>

#include "base/proteoform_reader.hpp"
#include "spec/peak.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_unexpected_change.hpp"
#include "prsmview/anno_expected_change.hpp"
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
      int file_num = getChildCount(root, "file");
      for (int i = 0; i < file_num; i++) {
        xercesc::DOMElement* file_element = getChildElement(root, "file", i);
        std::vector<std::string> file_info;
        file_info.push_back(prot::getChildValue(file_element,"xml",0));
        file_info.push_back(prot::getChildValue(file_element,"xsl",0));
        file_info.push_back(prot::getChildValue(file_element,"html",0));
        file_list.push_back(file_info);
      }
    }
    delete doc;
  }
  return file_list;
}

xercesc::DOMElement* genePrsmView(XmlDOMDocument* xml_doc,PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr){
  xercesc::DOMElement* element = xml_doc->createElement("prsm");
  std::string str = convertToString(prsm_ptr->getId());
  xml_doc->addElement(element, "prsm_id", str.c_str());
  if(prsm_ptr->getProbPtr().get()!=nullptr){
    str=convertToString(prsm_ptr->getProbPtr()->getPValue(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "p_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "p_value", "N/A");
  }
  if(prsm_ptr->getProbPtr().get()!=nullptr){
    str=convertToString(prsm_ptr->getProbPtr()->getEValue(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "e_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "e_value", "N/A");
  }
  double fdr = prsm_ptr->getFdr();
  if (fdr >= 0) {
    str=convertToString(prsm_ptr->getFdr(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "fdr", str.c_str());
  }
  else {
    xml_doc->addElement(element, "fdr", "N/A");
  }
  str=convertToString((int)prsm_ptr->getMatchFragNum());
  xml_doc->addElement(element, "matched_fragment_number", str.c_str());
  str=convertToString((int)prsm_ptr->getMatchPeakNum());
  xml_doc->addElement(element, "matched_peak_number", str.c_str());

  xercesc::DOMElement* ms_element = xml_doc->createElement("ms");
  xercesc::DOMElement* ms_header_element = xml_doc->createElement("ms_header");
  ms_element->appendChild(ms_header_element);
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  std::string spec_ids;
  std::string spec_scans;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + std::to_string(deconv_ms_ptr_vec[i]->getHeaderPtr()->getId()) + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getHeaderPtr()->getScansString() + " ";
  }
  boost::algorithm::trim(spec_ids);
  boost::algorithm::trim(spec_scans);
  xml_doc->addElement(ms_header_element, "ids", spec_ids.c_str());
  xml_doc->addElement(ms_header_element, "scans", spec_scans.c_str());
  int pos = 4;
  double precursor_mass = prsm_ptr->getOriPrecMass();
  str=convertToString(precursor_mass, pos);
  xml_doc->addElement(ms_header_element, "precursor_mono_mass", str.c_str());
  int precursor_charge = deconv_ms_ptr_vec[0]->getHeaderPtr()->getPrecCharge();
  str=convertToString(precursor_charge);
  xml_doc->addElement(ms_header_element, "precursor_charge", str.c_str());
  double precursor_mz = compMonoMz(precursor_mass, precursor_charge); 
  str=convertToString(precursor_mz, pos);
  xml_doc->addElement(ms_header_element, "precursor_mz", str.c_str());

  //peaks to view
  xercesc::DOMElement* peaks = xml_doc->createElement("peaks");
  ms_element->appendChild(peaks);
  ExtendMsPtrVec refine_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  for (size_t s = 0; s < deconv_ms_ptr_vec.size(); s++) {
    //get ion_pair
    PeakIonPairPtrVec pair_ptrs =  getPeakIonPairs (prsm_ptr->getProteoformPtr(), 
                                                    refine_ms_ptr_vec[s],
                                                    mng_ptr->min_mass_);
    //LOG_DEBUG("pair completed");
    for(size_t i=0;i< deconv_ms_ptr_vec[s]->size();i++){
      xercesc::DOMElement* peak_element = xml_doc->createElement("peak");
      peaks->appendChild(peak_element);
      str = convertToString(deconv_ms_ptr_vec[s]->getHeaderPtr()->getId());
      xml_doc->addElement(peak_element, "spec_id", str.c_str());
      DeconvPeakPtr peak_ptr = deconv_ms_ptr_vec[s]->getPeakPtr(i);
      str=convertToString(peak_ptr->getId());
      xml_doc->addElement(peak_element, "peak_id", str.c_str());
      double mass = peak_ptr->getPosition();
      int charge = peak_ptr->getCharge();
      str=convertToString(mass, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mass", str.c_str());
      double mz = compMonoMz(mass, charge);
      str=convertToString(mz, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mz", str.c_str());
      str=convertToString(peak_ptr->getIntensity(), mng_ptr->decimal_point_num_);
      xml_doc->addElement(peak_element, "intensity", str.c_str());
      str=convertToString(charge);
      xml_doc->addElement(peak_element, "charge", str.c_str());
      int spec_id = deconv_ms_ptr_vec[s]->getHeaderPtr()->getId(); 
      PeakIonPairPtrVec selected_pair_ptrs = getMatchedPairs(pair_ptrs, spec_id,  
                                                             peak_ptr->getId());
      if(selected_pair_ptrs.size()>0){
        int match_ions_number = selected_pair_ptrs.size();
        str=convertToString(match_ions_number);
        xml_doc->addElement(peak_element, "matched_ions_num", str.c_str());
        xercesc::DOMElement* mi_element = xml_doc->createElement("matched_ions");
        peak_element->appendChild(mi_element);
        for(size_t j=0;j< selected_pair_ptrs.size();j++){
          selected_pair_ptrs[j]->appendTheoPeakToXml(xml_doc,mi_element);
        }
      }
    }
  }
  element->appendChild(ms_element);
  //LOG_DEBUG("ms completed");

  //proteoform to view
  double err = prsm_ptr->getOriPrecMass() * 
      mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
  xercesc::DOMElement* prot_element = geneProteinView(xml_doc, prsm_ptr, mng_ptr, err);
  element->appendChild(prot_element);
  //LOG_DEBUG("protein view completed");

  return element;
  }


xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc,
                                     PrsmPtr prsm_ptr,
                                     PrsmViewMngPtr mng_ptr, double err) {
  xercesc::DOMElement* prot_element = xml_doc->createElement("annotated_protein");
  ProteoformPtr proteoform_ptr = prsm_ptr->getProteoformPtr();
  std::string str=convertToString(proteoform_ptr->getSeqId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str=convertToString(proteoform_ptr->getSpeciesId());
  xml_doc->addElement(prot_element, "proteoform_id", str.c_str());
  str=proteoform_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str=proteoform_ptr->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  double mass = proteoform_ptr->getMass();
  str=convertToString(mass, mng_ptr->decimal_point_num_);
  xml_doc->addElement(prot_element, "proteoform_mass", str.c_str());
  str=convertToString(proteoform_ptr->getProtModPtr()->getPtmPtr()->isAcetylation());
  xml_doc->addElement(prot_element, "n_acetylation", str.c_str());
  int unexpected_change_number = proteoform_ptr->getUnexpectedChangeNum();
  str=convertToString(unexpected_change_number);
  xml_doc->addElement(prot_element, "unexpected_change_number", str.c_str());

  ChangePtrVec change_ptrs = proteoform_ptr->getChangePtrVec(); 
  std::sort(change_ptrs.begin(),change_ptrs.end(),compareChangeTypeUpPosUp);

  xercesc::DOMElement* anno_element = xml_doc->createElement("annotation");
  prot_element->appendChild(anno_element);
  str=convertToString(proteoform_ptr->getDbResSeqPtr()->getLen());
  xml_doc->addElement(anno_element, "protein_length", str.c_str());

  str=convertToString(proteoform_ptr->getStartPos());
  xml_doc->addElement(anno_element, "first_residue_position", str.c_str());
  str=convertToString(proteoform_ptr->getEndPos());
  xml_doc->addElement(anno_element, "last_residue_position", str.c_str());
  str=convertToString(proteoform_ptr->getProtModPtr()->getPtmPtr()->isAcetylation());

  //LOG_DEBUG("summary completed");

  AnnoCleavagePtrVec cleavage_ptrs = getProteoCleavage(prsm_ptr, mng_ptr->min_mass_);

  //LOG_DEBUG("cleavage completed");

  int prot_len = proteoform_ptr->getDbResSeqPtr()->getLen();
  // obtain residue_ptrs 
  AnnoResiduePtrVec res_ptrs;
  for(int i=0;i< prot_len;i++){
    res_ptrs.push_back(AnnoResiduePtr(new AnnoResidue(proteoform_ptr->getDbResSeqPtr()->getResiduePtr(i), i)));
  }

  //LOG_DEBUG("residue completed");
  // add information for N-terminal truncation
  int start_pos = proteoform_ptr->getStartPos();
  for (int i =0; i < start_pos; i++) { 
    cleavage_ptrs[i]->setType(CLEAVAGE_TYPE_N_TRUNCATION);
    res_ptrs[i]->setType(ANNO_RESIDUE_TYPE_N_TRUNCATION);
  }

  if (start_pos > 0) {
    cleavage_ptrs[start_pos]->setType(CLEAVAGE_TYPE_SEQ_START);
  }
  //LOG_DEBUG("n-trunc completed");

  // add information for C-terminal truncation
  int end_pos = proteoform_ptr->getEndPos();
  if (end_pos < prot_len - 1) {
    cleavage_ptrs[end_pos + 1]->setType(CLEAVAGE_TYPE_SEQ_END);
  }

  for (int i = end_pos + 1; i < prot_len; i++) {
    cleavage_ptrs[i+1]->setType(CLEAVAGE_TYPE_C_TRUNCATION);
    res_ptrs[i]->setType(ANNO_RESIDUE_TYPE_C_TRUNCATION);
  }
  //LOG_DEBUG("c-trunc completed");

  AnnoUnexpectedChangePtrVec unexpected_change_ptrs;
  AnnoExpectedChangePtrVec expected_change_ptrs;
  int unexpected_shift_color = 0;
  int last_right = -1;
  for (size_t i = 0; i < change_ptrs.size(); i++) {
    // if the mass is less than 1 Da
    if (std::abs(change_ptrs[i]->getMassShift()) <= 1 + err)
        continue;
    // add information for known changes 
    int left_db_bp = change_ptrs[i]->getLeftBpPos() + start_pos;
    int right_db_bp = change_ptrs[i]->getRightBpPos() + start_pos;
    double shift = change_ptrs[i]->getMassShift();
    if (change_ptrs[i]->getChangeType() != UNEXPECTED_CHANGE) { 
      res_ptrs[left_db_bp]->setType(ANNO_RESIDUE_TYPE_KNOWN_CHANGE);

      AnnoExpectedChangePtr existing_ptr 
          = findExpectedChange(expected_change_ptrs, change_ptrs[i]->getChangeType(), change_ptrs[i]->getPtmPtr());
      if (existing_ptr == nullptr) {
        existing_ptr = AnnoExpectedChangePtr(new AnnoExpectedChange(change_ptrs[i]->getChangeType(), change_ptrs[i]->getPtmPtr()));
        expected_change_ptrs.push_back(existing_ptr);
      }
      std::string acid_letter = proteoform_ptr->getDbResSeqPtr()->getResiduePtr(left_db_bp)->getAcidPtr()->getOneLetter();
      existing_ptr->addOccurence(left_db_bp, acid_letter);
    }
    else {
      if (left_db_bp == right_db_bp) {
        int this_left = left_db_bp * 2;
        if (this_left > last_right + 1) {
          AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(last_right + 1, this_left - 1, 0, -1, "EMPTY"));
          unexpected_change_ptrs.push_back(anno_change_ptr);
        }
        int this_right = right_db_bp * 2;
        AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(this_left , this_right, shift, unexpected_shift_color, "SHIFT"));
        anno_change_ptr->setPtmPtr(change_ptrs[i]->getPtmPtr());
        std::string anno_info = "PTM: ";
        if (change_ptrs[i]->getPtmPtr() == nullptr) {
            anno_info += "Unknown";
        } else {
            anno_info += change_ptrs[i]->getPtmPtr()->getName();
        }

        for (int k = left_db_bp; k <= right_db_bp; k++) {
            std::string acid_letter = proteoform_ptr->getDbResSeqPtr()
                                      ->getResiduePtr(k)->getAcidPtr()->getOneLetter();
            anno_change_ptr->addOccurence(k, acid_letter);
            res_ptrs[k]->setPossiblePosColor(1);
            res_ptrs[k]->setAnno(anno_info);
        }


        unexpected_change_ptrs.push_back(anno_change_ptr);
        last_right = this_right;
        cleavage_ptrs[left_db_bp]->setUnexpectedChange(true);
        cleavage_ptrs[left_db_bp]->setUnexpectedChangeColor(unexpected_shift_color);;
      }
      else {
        int this_left = left_db_bp * 2 + 1;
        if (this_left > last_right + 1) {
          AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(last_right + 1, this_left - 1, 0, -1, "EMPTY"));
          unexpected_change_ptrs.push_back(anno_change_ptr);
        }
        int this_right = right_db_bp * 2 - 1;
        AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(this_left, this_right, shift, unexpected_shift_color, "SHIFT"));

        anno_change_ptr->setPtmPtr(change_ptrs[i]->getPtmPtr());
        std::string anno_info = "PTM: ";
        if (change_ptrs[i]->getPtmPtr() != nullptr) {
            anno_info += change_ptrs[i]->getPtmPtr()->getName() + "\n";
            std::vector<double> scr = change_ptrs[i]->getScr();
            for (int k = left_db_bp; k < right_db_bp; k++) {
                if (scr[k - left_db_bp] > 0) {
                    std::string acid_letter = proteoform_ptr->getDbResSeqPtr()
                        ->getResiduePtr(k)->getAcidPtr()->getOneLetter();
                    anno_info += "Site: " + acid_letter + std::to_string(k) + " ";
                    anno_info += "Confidence: "
                        + convertToString(scr[k - left_db_bp] * 100, 2) + "%\n";
                }
            }
            for (int k = left_db_bp; k < right_db_bp; k++) {
                if (scr[k - left_db_bp] > 0) {
                    std::string acid_letter = proteoform_ptr->getDbResSeqPtr()
                        ->getResiduePtr(k)->getAcidPtr()->getOneLetter();
                    anno_change_ptr->addOccurence(k, acid_letter);
                    res_ptrs[k]->setPossiblePosColor(1);
                    res_ptrs[k]->setAnno(anno_info);
                }
            }
        } else {
            anno_info += "Unknown\n";
            std::string acid_letter = proteoform_ptr->getDbResSeqPtr()
                ->getResiduePtr(left_db_bp)->getAcidPtr()->getOneLetter();
            anno_change_ptr->addOccurence(left_db_bp, acid_letter);
            anno_info += "Region: " + acid_letter + std::to_string(left_db_bp) + " - ";
            acid_letter = proteoform_ptr->getDbResSeqPtr()->getResiduePtr(
                    right_db_bp - 1)->getAcidPtr()->getOneLetter();
            anno_info += acid_letter + std::to_string(right_db_bp - 1);
            anno_change_ptr->addOccurence(right_db_bp - 1, acid_letter);
            double scr_sum = 0.0;
            std::vector<double> scr = change_ptrs[i]->getScr();
            for (int k = left_db_bp; k < right_db_bp; k++) {
                scr_sum += scr[k - left_db_bp];
            }
            anno_info += " Confindence: " + convertToString(scr_sum * 100, 2) + "%\n";
            for (int k = left_db_bp; k < right_db_bp; k++) {
                res_ptrs[k]->setPossiblePosColor(1);
                res_ptrs[k]->setAnno(anno_info);
            }

        }

        unexpected_change_ptrs.push_back(anno_change_ptr);
        last_right = this_right;

        for (int j = left_db_bp; j < right_db_bp - 1; j++) {
          res_ptrs[j]->setUnexpectedChange(true);
          res_ptrs[j]->setUnexpectedChangeColor(unexpected_shift_color);;
          cleavage_ptrs[j+1]->setUnexpectedChange(true);
          cleavage_ptrs[j+1]->setUnexpectedChangeColor(unexpected_shift_color);;
        }
        res_ptrs[right_db_bp-1]->setUnexpectedChange(true);
        res_ptrs[right_db_bp-1]->setUnexpectedChangeColor(unexpected_shift_color);;
      }
      unexpected_shift_color = (unexpected_shift_color) + 1 % 2;
    }
  }
  // last annochange
  AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(last_right + 1, std::numeric_limits<int>::max(), 0, -1, "EMPTY"));
  unexpected_change_ptrs.push_back(anno_change_ptr);
  /* remove EMPTY_CHANGES */
  AnnoUnexpectedChangePtr non_empty_ptr;
  AnnoUnexpectedChangePtrVec short_unexpected_change_ptrs;
  for (size_t i = 0; i < unexpected_change_ptrs.size(); i++) {
    AnnoUnexpectedChangePtr cur_change_ptr = unexpected_change_ptrs[i];
    if (cur_change_ptr->getChangeType() != "EMPTY") {
      short_unexpected_change_ptrs.push_back(cur_change_ptr);
      non_empty_ptr = cur_change_ptr;
    }
    else {
      if (non_empty_ptr == nullptr) {
        // first empty segment is kept 
        short_unexpected_change_ptrs.push_back(cur_change_ptr);
      }
      else {
        non_empty_ptr->setRightPos(cur_change_ptr->getRightPos());
      }
    }
  }

  //LOG_DEBUG("unexpected completed");

  for (size_t i = 0; i < res_ptrs.size(); i++) {
    res_ptrs[i]->appendViewXml(xml_doc, anno_element);
  }
  for (size_t i = 0; i < cleavage_ptrs.size(); i++) {
    cleavage_ptrs[i]->appendXml(xml_doc, anno_element);
  }
  for (size_t i = 0; i < short_unexpected_change_ptrs.size(); i++) {
    short_unexpected_change_ptrs[i]->appendXml(xml_doc, anno_element, mng_ptr->decimal_point_num_);
  }
  for (size_t i = 0; i < expected_change_ptrs.size(); i++) {
    expected_change_ptrs[i]->appendXml(xml_doc, anno_element);
  }
  return prot_element;
}



xercesc::DOMElement* proteoformToXml(XmlDOMDocument* xml_doc, const PrsmPtrVec &prsm_ptrs, 
                                     PrsmViewMngPtr mng_ptr){
  xercesc::DOMElement* proteoform_element = xml_doc->createElement("compatible_proteoform");
  std::string str=convertToString(prsm_ptrs[0]->getProteoformPtr()->getSeqId());
  xml_doc->addElement(proteoform_element, "sequence_id", str.c_str());
  str=prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(proteoform_element, "sequence_name", str.c_str());
  str=prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(proteoform_element, "sequence_description", str.c_str());
  str=convertToString(prsm_ptrs[0]->getProteoformPtr()->getSpeciesId());
  xml_doc->addElement(proteoform_element, "proteoform_id", str.c_str());
  int count = prsm_ptrs.size();
  str=convertToString(count);
  xml_doc->addElement(proteoform_element, "prsm_number", str.c_str());
  for(size_t i=0;i<prsm_ptrs.size();i++){
    proteoform_element->appendChild(genePrsmView(xml_doc,prsm_ptrs[i], mng_ptr));
  }
  return proteoform_element;
}

xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  ProteoformPtr proteo_ptr,
                                  const std::vector<int> &species_ids,
                                  PrsmViewMngPtr mng_ptr){
  xercesc::DOMElement* prot_element = xml_doc->createElement("protein");
  std::string str=convertToString(proteo_ptr->getSeqId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str=proteo_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str=proteo_ptr->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  int count = species_ids.size();
  str=convertToString(count);
  xml_doc->addElement(prot_element, "compatible_proteoform_number", str.c_str());
  for(size_t i=0;i<species_ids.size();i++){
    PrsmPtrVec select_prsm_ptrs = selectSpeciesPrsms(prsm_ptrs,species_ids[i]);
    std::sort(select_prsm_ptrs.begin(),select_prsm_ptrs.end(),prsmEValueUp);
    prot_element->appendChild(proteoformToXml(xml_doc,select_prsm_ptrs, mng_ptr));
  }
  return prot_element;
}

PrsmPtr getBestEValuePrsmPtr (ProteoformPtr proteo_ptr, const PrsmPtrVec &prsm_ptrs) {
  PrsmPtr best_ptr(nullptr);
  double best_evalue = std::numeric_limits<double>::max();
  int seq_id = proteo_ptr->getDbResSeqPtr()->getId();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getId() == seq_id && 
        prsm_ptrs[i]->getEValue() < best_evalue) {
      best_evalue = prsm_ptrs[i]->getEValue();
      best_ptr = prsm_ptrs[i];
    }
  }
  return best_ptr;
}

inline bool evalueCompare(const std::pair<ProteoformPtr, double> &a, const std::pair<ProteoformPtr, double> &b) {
    return a.second < b.second;
}


xercesc::DOMElement* allProteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  PrsmViewMngPtr mng_ptr){
  xercesc::DOMElement* prot_elements = xml_doc->createElement("proteins");
  // sort 
  ProteoformReader reader(mng_ptr->prsm_para_ptr_->getSearchDbFileName());
  ResiduePtrVec residue_ptr_vec = mng_ptr->prsm_para_ptr_->getFixModResiduePtrVec();
  ProteoformPtr proteo_ptr = reader.getNextProteoformPtr(residue_ptr_vec);
  std::vector<std::pair<ProteoformPtr, double>> proteo_evalues;

  while (proteo_ptr != nullptr) {
    PrsmPtr best_ptr = getBestEValuePrsmPtr (proteo_ptr, prsm_ptrs);
    if (best_ptr != nullptr) {
      std::pair<ProteoformPtr, double> cur_proteo_evalue(proteo_ptr, best_ptr->getEValue());
      proteo_evalues.push_back(cur_proteo_evalue);
    }
    proteo_ptr = reader.getNextProteoformPtr(residue_ptr_vec);
  }
  std::sort(proteo_evalues.begin(), proteo_evalues.end(), evalueCompare);
  
  for(size_t i=0;i<proteo_evalues.size();i++){
    std::vector<int> species_ids = getSpeciesIds(prsm_ptrs,proteo_evalues[i].first->getDbResSeqPtr()->getId());
    if(species_ids.size()>0){
      prot_elements->appendChild(proteinToXml(xml_doc,prsm_ptrs,proteo_evalues[i].first,species_ids, mng_ptr));
    }
  }
  return prot_elements;
}

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

}
