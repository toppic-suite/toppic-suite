/*
 * anno_view.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: xunlikun
 */

#include <xpp/anno_view.hpp>

namespace prot{

xercesc::DOMElement* AnnoView::geneFileList(XmlDOMDocument* xml_doc){
  xercesc::DOMElement* element = xml_doc->createElement("file_list");
  for(unsigned int i=0;i<file_list_.size();i++){
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

xercesc::DOMElement* genePrsmView(XmlDOMDocument* xml_doc,PrsmPtr prsm, double min_mass){
  int pos = 4;
  xercesc::DOMElement* element = xml_doc->createElement("prsm");
  std::string str = convertToString(prsm->getId());
  xml_doc->addElement(element, "prsm_id", str.c_str());
//  str = convertToString(prsm->getSpectrumId());
//  xml_doc->addElement(element, "spectrum_id", str.c_str());
//  xml_doc->addElement(element, "spectrum_scan", prsm->getSpectrumScan().c_str());
  if(prsm->getProbPtr().get()!=nullptr && prsm->getProbPtr()->getPValue() != -std::numeric_limits<double>::max()){
    str=convertToString(prsm->getProbPtr()->getPValue(),pos-2);
    xml_doc->addElement(element, "p_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "p_value", "-1");
  }
  if(prsm->getProbPtr().get()!=nullptr && prsm->getProbPtr()->getEValue() != -std::numeric_limits<double>::max()){
    str=convertToString(prsm->getProbPtr()->getEValue(),pos-2);
    xml_doc->addElement(element, "e_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "e_value", "-1");
  }
  str=convertToString(prsm->getFdr(),pos-2);
  xml_doc->addElement(element, "fdr", str.c_str());
  str=convertToString(prsm->getMatchFragNum(),pos-4);
  xml_doc->addElement(element, "matched_fragment_number", str.c_str());
  str=convertToString(prsm->getMatchPeakNum(),pos-4);
  xml_doc->addElement(element, "matched_peak_number", str.c_str());
//  str=convertToString(prsm->getOriPrecMass());
//  xml_doc->addElement(element, "precursor_mass", str.c_str());
  str=convertToString(prsm->getAdjustedPrecMass(),pos);
  xml_doc->addElement(element, "adjusted_precursor_mass", str.c_str());
  str=convertToString(prsm->getCalibration());
  xml_doc->addElement(element, "calibration", str.c_str());

  //get ion_pair
  PeakIonPairPtrVec pairs =  getPeakIonPairs (prsm->getProteoformPtr(), 
                                              prsm->getRefineMs(),
                                              min_mass);
  //peaks to view
  xercesc::DOMElement* ms_element = xml_doc->createElement("ms");
  prsm->getDeconvMsPtr()->getHeaderPtr()->appendXml(xml_doc,ms_element);//attention
  xercesc::DOMElement* peaks = xml_doc->createElement("peaks");
  ms_element->appendChild(peaks);
  for(unsigned int i=0;i<prsm->getDeconvMsPtr()->size();i++){
    xercesc::DOMElement* peak = xml_doc->createElement("peak");
    peaks->appendChild(peak);
    DeconvPeakPtr dp = prsm->getDeconvMsPtr()->getPeakPtr(i);
    str=convertToString(dp->getId());
    xml_doc->addElement(peak, "id", str.c_str());
    double mass = dp->getPosition();
    int charge = dp->getCharge();
    str=convertToString(mass,pos);
    xml_doc->addElement(peak, "monoisotopic_mass", str.c_str());
    str=convertToString(mass/charge+MassConstant::getProtonMass(),pos);
    xml_doc->addElement(peak, "monoisotopic_mz", str.c_str());
    str=convertToString(dp->getIntensity(),pos-2);
    xml_doc->addElement(peak, "intensity", str.c_str());
    str=convertToString(charge);
    xml_doc->addElement(peak, "charge", str.c_str());
    PeakIonPairPtrVec selected_pairs = getMatchedPairs(pairs,dp->getId());
    if(selected_pairs.size()>0){
      int match_ions_number = selected_pairs.size();
      str=convertToString(match_ions_number);
      xml_doc->addElement(peak, "matched_ions_num", str.c_str());
      xercesc::DOMElement* mi_element = xml_doc->createElement("matched_ions");
      peak->appendChild(mi_element);
      for(unsigned int j=0;j< selected_pairs.size();j++){
        selected_pairs[j]->appendIonToXml(xml_doc,mi_element);
      }
    }
  }
  element->appendChild(ms_element);

  //proteoform to view
  xercesc::DOMElement* prot_element = geneProteinView(xml_doc,
                                                      prsm->getProteoformPtr(),
                                                      prsm->getRefineMs(),
                                                      min_mass);
  element->appendChild(prot_element);

  return element;
}

//xercesc::DOMElement* genePrsmViewAS7(XmlDOMDocument* xml_doc,PrsmPtr prsm){
//  xercesc::DOMElement* element = xml_doc->createElement("prsm");
//  std::string str = convertToString(prsm->getId());
//  xml_doc->addElement(element, "prsm-id", str.c_str());
//
//  str=convertToString(prsm->getProteoformPtr()->getDbResSeqPtr()->getId());
//  xml_doc->addElement(element, "spectrum-id", str.c_str());
//  str=convertToString(prsm->getAdjustedPrecMass());
//  xml_doc->addElement(element, "spectrum-adjusted-mass", str.c_str());
//  str=prsm->getProteoformPtr()->getDbResSeqPtr()->getName();
//  xml_doc->addElement(element, "sequence-name", str.c_str());
//  str=convertToString(prsm->getProteoformPtr()->getSeqId());
//  xml_doc->addElement(element, "peptide-id", str.c_str());
//  //duplicate-score,used by peptides.xsl
//  //unique-score
//  if(prsm->getProbPtr()->getPValue() != -std::numeric_limits<double>::max()){
//    str=convertToString(prsm->getProbPtr()->getPValue());
//    xml_doc->addElement(element, "p_value", str.c_str());
//  }
//  else{
//    xml_doc->addElement(element, "p_value", "-1");
//  }
//  if(prsm->getProbPtr()->getEValue() != -std::numeric_limits<double>::max()){
//    str=convertToString(prsm->getProbPtr()->getEValue());
//    xml_doc->addElement(element, "e_value", str.c_str());
//  }
//  else{
//    xml_doc->addElement(element, "e_value", "-1");
//  }
//  str=convertToString(prsm->getFdr());
//  xml_doc->addElement(element, "fdr", str.c_str());
//  str=convertToString(prsm->getProteoformPtr()->getMass());
//  xml_doc->addElement(element, "sequence-mass", str.c_str());
//  //alignment type ,used by peptides.xsl
//
//  //tags
//  //sequence-name
//  //peakBreaks
//  str=convertToString(prsm->getMatchFragNum());
//  xml_doc->addElement(element, "matched_fragment_number", str.c_str());
//  str=convertToString(prsm->getMatchPeakNum());
//  xml_doc->addElement(element, "matched_peak_number", str.c_str());
//
//  str=convertToString(prsm->getCalibration());
//  xml_doc->addElement(element, "calibration", str.c_str());
//
//  //get ion_pair
//  PeakIonPairPtrVec pairs;
//  getPeakIonPairs (prsm->getProteoformPtr(), prsm->getRefineMs(),
//                   prsm->getMinMass(), pairs);
//  //peaks to view
//  xercesc::DOMElement* ms_element = xml_doc->createElement("ms");
//  prsm->getDeconvMsPtr()->getHeaderPtr()->appendXml(xml_doc,ms_element);//attention
//  xercesc::DOMElement* peaks = xml_doc->createElement("peaks");
//  ms_element->appendChild(peaks);
//  for(unsigned int i=0;i<prsm->getDeconvMsPtr()->size();i++){
//    xercesc::DOMElement* peak = xml_doc->createElement("peak");
//    peaks->appendChild(peak);
//    DeconvPeakPtr dp = prsm->getDeconvMsPtr()->getPeakPtr(i);
//    str=convertToString(dp->getId());
//    xml_doc->addElement(peak, "id", str.c_str());
//    double mass = dp->getPosition();
//    int charge = dp->getCharge();
//    str=convertToString(mass);
//    xml_doc->addElement(peak, "monoisotopic_mass", str.c_str());
//    str=convertToString(mass/charge+MassConstant::getProtonMass());
//    xml_doc->addElement(peak, "monoisotopic_mz", str.c_str());
//    str=convertToString(dp->getIntensity());
//    xml_doc->addElement(peak, "intensity", str.c_str());
//    str=convertToString(charge);
//    xml_doc->addElement(peak, "charge", str.c_str());
//    PeakIonPairPtrVec selected_pairs;
//    getMatchedPairs(pairs,dp->getId(),selected_pairs);
//    if(selected_pairs.size()>0){
//      xercesc::DOMElement* mi_element = xml_doc->createElement("matched_ions");
//      peak->appendChild(mi_element);
//      for(unsigned int j=0;j< selected_pairs.size();j++){
//        selected_pairs[j]->appendIonToXml(xml_doc,mi_element);
//      }
//    }
//  }
//  element->appendChild(ms_element);
//
//  //proteoform to view
//  xercesc::DOMElement* prot_element = geneProteinView(xml_doc,
//                                                      prsm->getProteoformPtr(),
//                                                      prsm->getRefineMs(),
//                                                      prsm->getMinMass());
//  element->appendChild(prot_element);
//
//  return element;
//}
//
//xercesc::DOMElement* geneProteinViewAS7(XmlDOMDocument* xml_doc,
//                                     ProteoformPtr proteoform_ptr,
//                                     ExtendMsPtr refine_ms_three,
//                                     double min_mass){
//  xercesc::DOMElement* prot_element = xml_doc->createElement("protein");
//  std::string str=convertToString(proteoform_ptr->getDbResSeqPtr()->getId());
//  str=proteoform_ptr->getDbResSeqPtr()->getName();
//  xml_doc->addElement(prot_element, "protein-name", str.c_str());
//
//  xml_doc->addElement(prot_element, "protein-id", str.c_str());
//  str=convertToString(proteoform_ptr->getSpeciesId());
//  xml_doc->addElement(prot_element, "species_id", str.c_str());
//
//  double mass = proteoform_ptr->getMass();
//  str=convertToString(mass);
//  xml_doc->addElement(prot_element, "protein_mass", str.c_str());
//  str=convertToString(proteoform_ptr->getStartPos());
//  xml_doc->addElement(prot_element, "first_residue_position", str.c_str());
//  str=convertToString(proteoform_ptr->getEndPos());
//  xml_doc->addElement(prot_element, "last_residue_position", str.c_str());
//  str=convertToString(proteoform_ptr->getProtModPtr()->getPtmPtr()->isAcetylation());
//  xml_doc->addElement(prot_element, "protein_mass", str.c_str());
//  int know_shift_number = proteoform_ptr->getChangePtrVec().size()-proteoform_ptr->getUnexpectedChangeNum();
//  str=convertToString(know_shift_number);
//  xml_doc->addElement(prot_element, "know_shift_number", str.c_str());
//  for(unsigned int i=0;i<proteoform_ptr->getChangePtrVec().size();i++){
//    proteoform_ptr->getChangePtrVec()[i]->appendXml(xml_doc,prot_element);//attention
//  }
//  xercesc::DOMElement* annotation_element = xml_doc->createElement("annotation");
//  prot_element->appendChild(annotation_element);
//  CleavagePtrVec cleavages = getProteoCleavage(proteoform_ptr,refine_ms_three,min_mass);
//  for(int i=0;i<proteoform_ptr->getResSeqPtr()->getLen();i++){
//    cleavages[i]->appendXml(xml_doc,annotation_element);
//    proteoform_ptr->getResSeqPtr()->getResiduePtr(i)->appendXml(xml_doc,annotation_element);
//  }
//  cleavages[cleavages.size()-1]->appendXml(xml_doc,annotation_element);
//  return prot_element;
//}

xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc,
                                     ProteoformPtr proteoform_ptr,
                                     ExtendMsPtr refine_ms_three,
                                     double min_mass) {
  int pos=4;
  xercesc::DOMElement* prot_element = xml_doc->createElement("annotated_protein");
  std::string str=convertToString(proteoform_ptr->getSeqId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str=convertToString(proteoform_ptr->getSpeciesId());
  xml_doc->addElement(prot_element, "species_id", str.c_str());
  str=proteoform_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  double mass = proteoform_ptr->getMass();
  str=convertToString(mass,pos);
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
//  for(unsigned int i=0;i<proteoform_ptr->getChangePtrVec().size();i++){
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
  for(unsigned int i=0;i<proteoform_ptr->getChangePtrVec().size();i++){
    if(proteoform_ptr->getChangePtrVec()[i]->getChangeType()!=UNEXPECTED_CHANGE){
      std::string abb_name = proteoform_ptr->getChangePtrVec()[i]->getPtmPtr()->getAbbrName();
      if(style[abb_name].compare("")==0){
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
  CleavagePtrVec cleavages = getProteoCleavage(proteoform_ptr,refine_ms_three,min_mass);
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

    if(i>0 && cleavages[i]->getType().compare("n_truncation")!=0 && cleavages[i-1]->getType().compare("n_truncation")==0){
      cleavages[i]->setTrunc("]");
    }

//    if(i>0)
//    std::cout<<(cleavages[i]->getType().compare("c_truncation")==0)<<(cleavages[i-1]->getType().compare("c_truncation")!=0)<<std::endl;
    if(i>0 && cleavages[i]->getType().compare("c_truncation")==0 && cleavages[i-1]->getType().compare("c_truncation")!=0){
      cleavages[i]->setTrunc("[");
    }
//    std::cout<<cleavages[i]->getTrunc()<<std::endl;

    ChangePtrVec change_list = proteoform_ptr->getChangePtrVec();
    for(unsigned int j=0;j<change_list.size();j++){
      if(change_list[j]->getLeftBpPos()+proteoform_ptr->getStartPos()-1<i&& i<change_list[j]->getRightBpPos()+proteoform_ptr->getStartPos())
      {
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
      if(change_list[j]->getLeftBpPos()==change_list[j]->getRightBpPos() and change_list[j]->getLeftBpPos()==i){
        cleavages[i]->setType("unexpected_shift");
        cleavages[i]->setShift(change_list[j]->getMassShift());
        if (j > 0 && change_list[j]->getLeftBpPos() - change_list[j - 1]->getLeftBpPos() <= 5) {
          display = 1 - display;
        } else {
          display = 0;
        }
        cleavages[i]->setDisplayPos(display);
      }
//      display_bg = 1-display_bg;
    }
    cleavages[i]->appendXml(xml_doc,annotation_element);
    cur_res->appendViewXml(xml_doc,annotation_element);
  }
  cleavages[cleavages.size()-1]->appendXml(xml_doc,annotation_element);
  return prot_element;
}

std::vector<int> getSpeciesIds(PrsmPtrVec prsms,int seq_id){
  std::vector<int> species;
  for(unsigned int i=0;i<prsms.size();i++){
    int new_id = prsms[i]->getProteoformPtr()->getSpeciesId();
    if(prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getId() == seq_id){
      bool flag= false;
      for(unsigned int j=0;j<species.size();j++){
        if(species[j]==new_id){
          flag=true;
        }
      }
      if(!flag){
        species.push_back(new_id);
      }
    }
  }
  return species;
}

std::vector<int> getSpeciesIds(PrsmPtrVec prsms){
//  std::cout<<prsms.size()<<std::endl;
  std::vector<int> species;
  for(unsigned int i=0;i<prsms.size();i++){
    bool find = false;
    for(unsigned int j=0;j<species.size();j++){
      if(species[j]==prsms[i]->getProteoformPtr()->getSpeciesId()){
        find = true;
      }
    }
    if(!find){
      species.push_back(prsms[i]->getProteoformPtr()->getSpeciesId());
    }
  }

  return species;
}

PrsmPtrVec selectSpeciesPrsms(PrsmPtrVec prsms,int species_id){
  PrsmPtrVec select_prsms;
  for(unsigned int i=0;i<prsms.size();i++){
    if(species_id == prsms[i]->getProteoformPtr()->getSpeciesId()){
      select_prsms.push_back(prsms[i]);
    }
  }
  return select_prsms;
}

xercesc::DOMElement* speciesToXml(XmlDOMDocument* xml_doc,PrsmPtrVec prsms, double min_mass){
  xercesc::DOMElement* species_element = xml_doc->createElement("species");
  std::string str=convertToString(prsms[0]->getProteoformPtr()->getSeqId());
  xml_doc->addElement(species_element, "sequence_id", str.c_str());
  str=prsms[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(species_element, "sequence_name", str.c_str());
  str=convertToString(prsms[0]->getProteoformPtr()->getSpeciesId());
  xml_doc->addElement(species_element, "species_id", str.c_str());
  int count = prsms.size();
  str=convertToString(count);
  xml_doc->addElement(species_element, "prsm_number", str.c_str());
  for(unsigned int i=0;i<prsms.size();i++){
    species_element->appendChild(genePrsmView(xml_doc,prsms[i], min_mass));
  }
  return species_element;
}

xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  PrsmPtrVec prsms,
                                  ProteoformPtr protein,
                                  std::vector<int> species,
                                  double min_mass){
  xercesc::DOMElement* prot_element = xml_doc->createElement("protein");
  std::string str=convertToString(protein->getSeqId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str=protein->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  int count = species.size();
  str=convertToString(count);
  xml_doc->addElement(prot_element, "species_number", str.c_str());
  for(unsigned int i=0;i<species.size();i++){
    PrsmPtrVec select_prsms = selectSpeciesPrsms(prsms,species[i]);
    std::sort(select_prsms.begin(),select_prsms.end(),prsmEValueUp);
    prot_element->appendChild(speciesToXml(xml_doc,select_prsms, min_mass));
  }
  return prot_element;
}

xercesc::DOMElement* allProteinToXml(XmlDOMDocument* xml_doc,
                                  PrsmPtrVec prsms,
                                  ProteoformPtrVec proteins,
                                  double min_mass){
  xercesc::DOMElement* prot_elements = xml_doc->createElement("proteins");
  for(unsigned int i=0;i<proteins.size();i++){
    std::vector<int> species = getSpeciesIds(prsms,proteins[i]->getDbResSeqPtr()->getId());
    if(species.size()>0){
      prot_elements->appendChild(proteinToXml(xml_doc,prsms,proteins[i],species, min_mass));
    }
  }
//  std::string str=convertToString(protein->getSeqId());
//  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
//  str=protein->getSeqName();
//  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
//  int count = species.size();
//  str=convertToString(count);
//  xml_doc->addElement(prot_element, "species_number", str.c_str());
//  for(unsigned int i=0;i<species.size();i++){
//    PrsmPtrVec select_prsms = selectSpeciesPrsms(prsms,species[i]);
//    std::sort(select_prsms.begin(),select_prsms.end(),prsmEValueDown);
//    prot_element->appendChild(speciesToXml(xml_doc,select_prsms));
//  }
  return prot_elements;
}

}
