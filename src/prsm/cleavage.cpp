#include "prsm/cleavage.hpp"

namespace prot {

Cleavage::Cleavage(int pos){
  pos_= pos;
  exist_n_ion_ = false;
  exist_c_ion_ = false;
  is_unexpected_change_ = false;
  unexpected_change_color_ = 0;
  type_ = CLEAVAGE_TYPE_NORMAL;
}

void Cleavage::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("cleavage");
  std::string str = convertToString(pos_);
  xml_doc->addElement(element, "position", str.c_str());
  str = type_;
  xml_doc->addElement(element, "cleavage_type", str.c_str());
  str = convertToString(is_unexpected_change_);
  xml_doc->addElement(element, "is_unexpected_change", str.c_str());
  str = convertToString(unexpected_change_color_);
  xml_doc->addElement(element, "unexpected_change_color", str.c_str());
  str = convertToString(exist_n_ion_);
  xml_doc->addElement(element, "exist_n_ion", str.c_str());
  str = convertToString(exist_c_ion_);
  xml_doc->addElement(element, "exist_c_ion", str.c_str());
  xercesc::DOMElement* peaks = xml_doc->createElement("matched_peaks");
  for(size_t i=0;i<pairs_.size();i++){
    pairs_[i]->appendRealPeakToXml(xml_doc,peaks);
  }
  element->appendChild(peaks);
  parent->appendChild(element);
}

CleavagePtrVec getProteoCleavage(ProteoformPtr prot_ptr,
                                 ExtendMsPtr ms_three_ptr,
                                 double min_mass){
  CleavagePtrVec cleavages;
  PeakIonPairPtrVec pairs = getPeakIonPairs (prot_ptr, ms_three_ptr, min_mass);

  PeakIonPairPtrVec2D peak_list;
  int prot_len = prot_ptr->getDbResSeqPtr()->getLen();
  std::vector<bool> n_ion (prot_len + 1, false);
  std::vector<bool> c_ion (prot_len + 1, false);
  for(int i=0; i< prot_len + 1; i++){
    PeakIonPairPtrVec temp;
    peak_list.push_back(temp);
  }
  for(size_t i=0;i<pairs.size();i++){
    int pos = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos()+prot_ptr->getStartPos();
    //LOG_DEBUG("start pos " << prot_ptr->getStartPos() << " pos " << pos);
    peak_list[pos].push_back(pairs[i]);
    if(pairs[i]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->isNTerm()){
      n_ion[pos] = true;
    }
    else{
      c_ion[pos] = true;
    }
  }

  for(int i=0;i< prot_len+1;i++){
    CleavagePtr cleavage = CleavagePtr(new Cleavage(i));
    cleavage->setPairs(peak_list[i]);
    cleavage->setExistNIon(n_ion[i]);
    cleavage->setExistCIon(c_ion[i]);
    cleavages.push_back(cleavage);
    //LOG_DEBUG("i  " << i << " n ion " << n_ion[i] << " c ion " << c_ion[i]);
  }
  return cleavages;
}


} /* namespace prot */
