
#include "spec/peak.hpp"
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
  //  str = convertToString(prsm->getSpectrumId());
  //  xml_doc->addElement(element, "spectrum_id", str.c_str());
  //  xml_doc->addElement(element, "spectrum_scan", prsm->getSpectrumScan().c_str());
  if(prsm_ptr->getProbPtr().get()!=nullptr){
    str=convertToString(prsm_ptr->getProbPtr()->getPValue(), mng_ptr->decimal_point_num);
    xml_doc->addElement(element, "p_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "p_value", "N/A");
  }
  if(prsm_ptr->getProbPtr().get()!=nullptr){
    str=convertToString(prsm_ptr->getProbPtr()->getEValue(), mng_ptr->decimal_point_num);
    xml_doc->addElement(element, "e_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "e_value", "N/A");
  }
  double fdr = prsm_ptr->getFdr();
  if (fdr < -1) {
    str=convertToString(prsm_ptr->getFdr(), mng_ptr->decimal_point_num);
    xml_doc->addElement(element, "fdr", str.c_str());
  }
  else {
    xml_doc->addElement(element, "fdr", "N/A");
  }
  str=convertToString((int)prsm_ptr->getMatchFragNum());
  xml_doc->addElement(element, "matched_fragment_number", str.c_str());
  str=convertToString((int)prsm_ptr->getMatchPeakNum());
  xml_doc->addElement(element, "matched_peak_number", str.c_str());
  //  str=convertToString(prsm->getOriPrecMass());
  //  xml_doc->addElement(element, "precursor_mass", str.c_str());
  str=convertToString(prsm_ptr->getAdjustedPrecMass(), mng_ptr->precise_point_num);
  xml_doc->addElement(element, "adjusted_precursor_mass", str.c_str());
  str=convertToString(prsm_ptr->getCalibration(), mng_ptr->precise_point_num);
  xml_doc->addElement(element, "calibration", str.c_str());

  //get ion_pair
  PeakIonPairPtrVec pair_ptrs =  getPeakIonPairs (prsm_ptr->getProteoformPtr(), 
                                                  prsm_ptr->getRefineMs(),
                                              mng_ptr->min_mass);
  //peaks to view
  xercesc::DOMElement* ms_element = xml_doc->createElement("ms");
  prsm_ptr->getDeconvMsPtr()->getHeaderPtr()->appendXml(xml_doc,ms_element);//attention
  xercesc::DOMElement* peaks = xml_doc->createElement("peaks");
  ms_element->appendChild(peaks);
  for(size_t i=0;i<prsm_ptr->getDeconvMsPtr()->size();i++){
    xercesc::DOMElement* peak_element = xml_doc->createElement("peak");
    peaks->appendChild(peak_element);
    DeconvPeakPtr peak_ptr = prsm_ptr->getDeconvMsPtr()->getPeakPtr(i);
    str=convertToString(peak_ptr->getId());
    xml_doc->addElement(peak_element, "id", str.c_str());
    double mass = peak_ptr->getPosition();
    int charge = peak_ptr->getCharge();
    str=convertToString(mass, mng_ptr->precise_point_num);
    xml_doc->addElement(peak_element, "monoisotopic_mass", str.c_str());
    double mz = compMonoMz(mass, charge);
    str=convertToString(mz, mng_ptr->precise_point_num);
    xml_doc->addElement(peak_element, "monoisotopic_mz", str.c_str());
    str=convertToString(peak_ptr->getIntensity(), mng_ptr->decimal_point_num);
    xml_doc->addElement(peak_element, "intensity", str.c_str());
    str=convertToString(charge);
    xml_doc->addElement(peak_element, "charge", str.c_str());
    PeakIonPairPtrVec selected_pair_ptrs = getMatchedPairs(pair_ptrs,peak_ptr->getId());
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
  element->appendChild(ms_element);

  //proteoform to view
  xercesc::DOMElement* prot_element = geneProteinView(xml_doc,
                                                      prsm_ptr->getProteoformPtr(),
                                                      prsm_ptr->getRefineMs(),
                                                      mng_ptr);
  element->appendChild(prot_element);

  return element;
  }


xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc,
                                     ProteoformPtr proteoform_ptr,
                                     ExtendMsPtr ms_three_ptr,
                                     PrsmViewMngPtr mng_ptr) {
  xercesc::DOMElement* prot_element = xml_doc->createElement("annotated_protein");
  std::string str=convertToString(proteoform_ptr->getSeqId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str=convertToString(proteoform_ptr->getSpeciesId());
  xml_doc->addElement(prot_element, "species_id", str.c_str());
  str=proteoform_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  double mass = proteoform_ptr->getMass();
  str=convertToString(mass, mng_ptr->decimal_point_num);
  xml_doc->addElement(prot_element, "protein_mass", str.c_str());
  str=convertToString(proteoform_ptr->getStartPos());
  xml_doc->addElement(prot_element, "first_residue_position", str.c_str());
  str=convertToString(proteoform_ptr->getEndPos());
  xml_doc->addElement(prot_element, "last_residue_position", str.c_str());
  str=convertToString(proteoform_ptr->getProtModPtr()->getPtmPtr()->isAcetylation());
  xml_doc->addElement(prot_element, "n_acetylation", str.c_str());
  int know_shift_number = proteoform_ptr->getChangePtrVec().size()-proteoform_ptr->getUnexpectedChangeNum();
  str=convertToString(know_shift_number);
  xml_doc->addElement(prot_element, "know_shift_number", str.c_str());
  str=convertToString(proteoform_ptr->getDbResSeqPtr()->getLen());
  xml_doc->addElement(prot_element, "db_acid_number", str.c_str());
//  for(size_t i=0;i<proteoform_ptr->getChangePtrVec().size();i++){
//    proteoform_ptr->getChangePtrVec()[i]->appendViewXml(xml_doc,prot_element);//attention
//  }
  xercesc::DOMElement* shift_list = xml_doc->createElement("shift_list");
  prot_element->appendChild(shift_list);
  std::vector<std::string> colors;
  colors.push_back("red");
  colors.push_back("blue");
  colors.push_back("green");
  colors.push_back("brown");
  colors.push_back("orange");
  std::map<std::string,std::string> style;
  int m=0;
  for(size_t i=0;i<proteoform_ptr->getChangePtrVec().size();i++){
    if(proteoform_ptr->getChangePtrVec()[i]->getChangeType()!=UNEXPECTED_CHANGE){
      std::string abb_name = proteoform_ptr->getChangePtrVec()[i]->getPtmPtr()->getAbbrName();
      if(style[abb_name] == ""){
        style[abb_name]=colors[m%colors.size()];
        xercesc::DOMElement* shift_element = xml_doc->createElement("shift");
        shift_list->appendChild(shift_element);
        xml_doc->addElement(shift_element, "type", abb_name.c_str());
        xml_doc->addElement(shift_element, "color", style[abb_name].c_str());
        str=convertToString(proteoform_ptr->getChangePtrVec()[i]->getChangeType());
        xml_doc->addElement(shift_element, "known_type", str.c_str());
      }
      m++;
    }
  }
  xercesc::DOMElement* annotation_element = xml_doc->createElement("annotation");
  prot_element->appendChild(annotation_element);
  CleavagePtrVec cleavages = getProteoCleavage(proteoform_ptr,ms_three_ptr,mng_ptr->min_mass);
  int display =0;
  //int display_bg =0;
  for(int i=0;i<proteoform_ptr->getDbResSeqPtr()->getLen();i++){
    cleavages[i]->setType("species");
    cleavages[i]->setTrunc("");
    AnnoResiduePtr cur_res = AnnoResiduePtr(new AnnoResidue(proteoform_ptr->getDbResSeqPtr()->getResiduePtr(i)));
    cur_res->setPos(i);
    cur_res->setType("normal");
    cur_res->setDisplayPos(0);
    cur_res->setIsModifyed(false);
    cur_res->setShift(0);
    cur_res->setShiftStyle("black");
    if(i<proteoform_ptr->getStartPos()){
      cur_res->setType("n_trunc");
      cleavages[i]->setType("n_truncation");
    }
    if(i>proteoform_ptr->getEndPos()){
      cur_res->setType("c_trunc");
      cleavages[i]->setType("c_truncation");
    }

    if(i>0 && cleavages[i]->getType() != "n_truncation" && cleavages[i-1]->getType() == "n_truncation") {
      cleavages[i]->setTrunc("]");
    } 

    if(i>0 && cleavages[i]->getType() == "c_truncation" && cleavages[i-1]->getType() != "c_truncation") {
      cleavages[i]->setTrunc("[");
    }

    ChangePtrVec change_list = proteoform_ptr->getChangePtrVec();
    for(size_t j=0;j<change_list.size();j++){
      if(change_list[j]->getLeftBpPos()+proteoform_ptr->getStartPos()-1<i && 
         i<change_list[j]->getRightBpPos()+proteoform_ptr->getStartPos()) {
        if(change_list[j]->getChangeType()==UNEXPECTED_CHANGE){
          cur_res->setType("unexpected_shift");
          cleavages[i]->setType("unexpected_shift");
          cur_res->setDisplayBg(j%2);
        }
        else{
          cur_res->setExpected(true);
          cleavages[i]->setType("expected_shift");
        }
        if(i==change_list[j]->getLeftBpPos()+proteoform_ptr->getStartPos()){
          if(change_list[j]->getChangeType()==UNEXPECTED_CHANGE){
            cur_res->setIsModifyed(true);
          }
          else{
            cur_res->setIsModifyed(false);
            std::string abb_name = change_list[j]->getPtmPtr()->getAbbrName();
            cur_res->setShiftStyle(style[abb_name]);
          }
          cur_res->setShift(change_list[j]->getMassShift());
          if(j>0&&change_list[j]->getLeftBpPos()-change_list[j-1]->getLeftBpPos()<=5){
            display=1-display;
          }
          else{
            display=0;
          }
          cur_res->setDisplayPos(display);
        }
      }
      if(change_list[j]->getLeftBpPos()==change_list[j]->getRightBpPos() 
         && change_list[j]->getLeftBpPos()==i){
        cleavages[i]->setType("unexpected_shift");
        cleavages[i]->setShift(change_list[j]->getMassShift());
        if (j > 0 && change_list[j]->getLeftBpPos() - change_list[j - 1]->getLeftBpPos() <= 5) {
          display = 1 - display;
        } else {
          display = 0;
        }
        cleavages[i]->setDisplayPos(display);
      }
    }
    cleavages[i]->appendXml(xml_doc,annotation_element);
    cur_res->appendViewXml(xml_doc,annotation_element);
  }
  cleavages[cleavages.size()-1]->appendXml(xml_doc,annotation_element);
  return prot_element;
}



xercesc::DOMElement* speciesToXml(XmlDOMDocument* xml_doc, const PrsmPtrVec &prsm_ptrs, 
                                  PrsmViewMngPtr mng_ptr){
  xercesc::DOMElement* species_element = xml_doc->createElement("species");
  std::string str=convertToString(prsm_ptrs[0]->getProteoformPtr()->getSeqId());
  xml_doc->addElement(species_element, "sequence_id", str.c_str());
  str=prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(species_element, "sequence_name", str.c_str());
  str=convertToString(prsm_ptrs[0]->getProteoformPtr()->getSpeciesId());
  xml_doc->addElement(species_element, "species_id", str.c_str());
  int count = prsm_ptrs.size();
  str=convertToString(count);
  xml_doc->addElement(species_element, "prsm_number", str.c_str());
  for(size_t i=0;i<prsm_ptrs.size();i++){
    species_element->appendChild(genePrsmView(xml_doc,prsm_ptrs[i], mng_ptr));
  }
  return species_element;
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
  int count = species_ids.size();
  str=convertToString(count);
  xml_doc->addElement(prot_element, "species_number", str.c_str());
  for(size_t i=0;i<species_ids.size();i++){
    PrsmPtrVec select_prsm_ptrs = selectSpeciesPrsms(prsm_ptrs,species_ids[i]);
    std::sort(select_prsm_ptrs.begin(),select_prsm_ptrs.end(),prsmEValueUp);
    prot_element->appendChild(speciesToXml(xml_doc,select_prsm_ptrs, mng_ptr));
  }
  return prot_element;
}

xercesc::DOMElement* allProteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  const ProteoformPtrVec &proteo_ptrs,
                                  PrsmViewMngPtr mng_ptr){
  xercesc::DOMElement* prot_elements = xml_doc->createElement("proteins");
  for(size_t i=0;i<proteo_ptrs.size();i++){
    std::vector<int> species_ids = getSpeciesIds(prsm_ptrs,proteo_ptrs[i]->getDbResSeqPtr()->getId());
    if(species_ids.size()>0){
      prot_elements->appendChild(proteinToXml(xml_doc,prsm_ptrs,proteo_ptrs[i],species_ids, mng_ptr));
    }
  }
//  std::string str=convertToString(protein->getSeqId());
//  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
//  str=protein->getSeqName();
//  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
//  int count = species.size();
//  str=convertToString(count);
//  xml_doc->addElement(prot_element, "species_number", str.c_str());
//  for(size_t i=0;i<species.size();i++){
//    PrsmPtrVec select_prsms = selectSpeciesPrsms(prsms,species[i]);
//    std::sort(select_prsms.begin(),select_prsms.end(),prsmEValueDown);
//    prot_element->appendChild(speciesToXml(xml_doc,select_prsms));
//  }
  return prot_elements;
}

}
