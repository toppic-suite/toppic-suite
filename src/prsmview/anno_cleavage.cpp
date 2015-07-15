#include "prsmview/anno_cleavage.hpp"

namespace prot {

AnnoCleavage::AnnoCleavage(int pos, const PeakIonPairPtrVec &pairs, 
                           bool exist_n_ion, bool exist_c_ion){
  pos_= pos;
  pairs_ = pairs;
  exist_n_ion_ = exist_n_ion;
  exist_c_ion_ = exist_c_ion;
  is_unexpected_change_ = false;
  unexpected_change_color_ = 0;
  type_ = CLEAVAGE_TYPE_NORMAL;
}

void AnnoCleavage::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
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

AnnoCleavagePtrVec getProteoCleavage(PrsmPtr prsm_ptr, double min_mass){
  AnnoCleavagePtrVec cleavages;
  ProteoformPtr proteo_ptr = prsm_ptr->getProteoformPtr();
  ExtendMsPtrVec refine_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  int prot_len = proteo_ptr->getDbResSeqPtr()->getLen();
  std::vector<bool> n_ion (prot_len + 1, false);
  std::vector<bool> c_ion (prot_len + 1, false);
  PeakIonPairPtrVec2D peak_list(prot_len + 1, PeakIonPairPtrVec(0));

  for (size_t m = 0; m < refine_ms_ptr_vec.size(); m++) {
    PeakIonPairPtrVec pairs = getPeakIonPairs (proteo_ptr, refine_ms_ptr_vec[m], 
                                               min_mass);
    for(size_t i=0; i<pairs.size(); i++){
      int pos = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos()+ proteo_ptr->getStartPos();
      //LOG_DEBUG("start pos " << prot_ptr->getStartPos() << " pos " << pos);
      peak_list[pos].push_back(pairs[i]);
      if(pairs[i]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->isNTerm()){
        n_ion[pos] = true;
      }
      else{
        c_ion[pos] = true;
      }
    }
  }

  for(int i=0;i< prot_len+1;i++){
    AnnoCleavagePtr cleavage = AnnoCleavagePtr(
        new AnnoCleavage(i, peak_list[i], n_ion[i], c_ion[i]));
    cleavages.push_back(cleavage);
    //LOG_DEBUG("i  " << i << " n ion " << n_ion[i] << " c ion " << c_ion[i]);
  }
  return cleavages;
}


} /* namespace prot */
