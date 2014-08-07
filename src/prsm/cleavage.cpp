#include "prsm/cleavage.hpp"

namespace prot {

Cleavage::Cleavage(int pos){
  pos_=pos;
  exist_n_ion_ = false;
  exist_c_ion_ = false;
  shift_= 0.0;
  display_pos_ = 0;
}

void Cleavage::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("character");
      std::string str = "cleavage";
      xml_doc->addElement(element, "type", str.c_str());
      str = type_;
      xml_doc->addElement(element, "cleavage_type", str.c_str());
      str = trunc_;
      xml_doc->addElement(element, "cleavage_trunc", str.c_str());
      str = convertToString(pos_);
      xml_doc->addElement(element, "position", str.c_str());
      str = convertToString(display_pos_);
      xml_doc->addElement(element, "display_position", str.c_str());
      str = convertToString(exist_n_ion_);
      xml_doc->addElement(element, "exist_n_ion", str.c_str());
      str = convertToString(exist_c_ion_);
      xml_doc->addElement(element, "exist_c_ion", str.c_str());
      str = convertToString(shift_,2);
      xml_doc->addElement(element, "shift_no_letter", str.c_str());
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
  std::vector<bool> n_ion;
  std::vector<bool> c_ion;
  for(size_t i=0; i<prot_ptr->getDbResSeqPtr()->getResidues().size()+1; i++){
    PeakIonPairPtrVec temp;
    peak_list.push_back(temp);
    n_ion.push_back(false);
    c_ion.push_back(false);
  }

  for(size_t i=0;i<pairs.size();i++){
    int pos = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos()+prot_ptr->getStartPos();
    peak_list[pos].push_back(pairs[i]);
    if(pairs[i]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->isNTerm()){
      n_ion[pos] = true;
    }
    else{
      c_ion[pos] = true;
    }
  }

  for(size_t i=0;i< prot_ptr->getDbResSeqPtr()->getResidues().size()+1;i++){
    CleavagePtr cleavage = CleavagePtr(new Cleavage(i));
    cleavage->setPairs(peak_list[i]);
    cleavage->setExistCIon(c_ion[i]);
    cleavage->setExistNIon(n_ion[i]);
    cleavages.push_back(cleavage);
  }
  return cleavages;
}


} /* namespace prot */
