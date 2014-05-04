/*
 * extend_peak.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include "spec/extend_peak.hpp"

namespace prot {

ExtendPeak::ExtendPeak(const DeconvPeakPtr &base_peak_ptr, 
                       double mono_mass,double score)
    :Peak(mono_mass,1.0){
  base_peak_ptr_ = base_peak_ptr;
  mono_mass_ = mono_mass;
  score_ = score;
  orig_tolerance_ = 0.0;
  reverse_tolerance_ = 0.0;
}

void ExtendPeak::appendXml(XmlDOMDocument* xml_doc,
                           xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("extend_peak");
  std::string str = convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = convertToString(mono_mass_);
  xml_doc->addElement(element, "mono_mass", str.c_str());
  str = convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  str = convertToString(orig_tolerance_);
  xml_doc->addElement(element, "orig_tolerance", str.c_str());
  str = convertToString(reverse_tolerance_);
  xml_doc->addElement(element, "reverse_tolerance_", str.c_str());
  base_peak_ptr_->appendXml(xml_doc,element);
  parent->appendChild(element);
}

//must deleted the ms after finished using
//ExtendMsPtr getMsThree(DeconvMsPtr deconv_ms,double delta,SpParaPtr sp_para);
ExtendMsPtr getMsThree(const DeconvMsPtr &deconv_ms, double delta,
                       const SpParaPtr &sp_para_ptr) {
  MsHeaderPtr header = getDeltaHeaderPtr(deconv_ms, delta);

  //private function getSpThreeExtendPeak in factory
  ExtendPeakPtrVec list;
  double ext_min_mass = sp_para_ptr->getExtendMinMass();
  std::vector<double> ext_offsets = sp_para_ptr->getExtendOffsets();
  for(unsigned int i =0; i<deconv_ms->size();i++){
    DeconvPeakPtr peak_ptr = deconv_ms->getPeakPtr(i);
    if(peak_ptr->getMonoMass() <= ext_min_mass) {
      double orig_mass = peak_ptr->getMonoMass();
      ExtendPeakPtr new_peak 
          = ExtendPeakPtr(new ExtendPeak(peak_ptr,orig_mass,1.0));
      list.push_back(new_peak);
    }
    else{
      for(unsigned int j=0;j < ext_offsets.size();j++){
        double mass = peak_ptr->getMonoMass() + ext_offsets[j];
        ExtendPeakPtr new_peak 
            = ExtendPeakPtr(new ExtendPeak(peak_ptr, mass,1.0));
        list.push_back(new_peak);
      }
    }
  }
  //filter extend_peak
  ExtendPeakPtrVec list_filtered;
  double min_mass = sp_para_ptr->getMinMass();
  double prec_mono_mass = header->getPrecMonoMass();
  for(unsigned int i =0; i < list.size();i++){
    double mass = list[i]->getPosition();
    if(mass >= min_mass && mass <= prec_mono_mass - min_mass){
      list_filtered.push_back(list[i]);
    }
  }
  //end filter extend_peak function and result = list_filtered

  std::sort(list_filtered.begin(),list_filtered.end(),extendPeakUp);
  //end getSpThreeExtendPeak and result = list_filtered;

  //function msThreeSetTolerance
  PeakTolerancePtr peak_tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  for (unsigned int i = 0; i < list_filtered.size();i++){
    double mass = list_filtered[i]->getBasePeakPtr()->getMonoMass();
    double ori_tole = peak_tole_ptr->compStrictErrorTole(mass);
    list_filtered[i]->setOrigTolerance(ori_tole);
    double reve_tole 
        = peak_tole_ptr->compRelaxErrorTole(mass, prec_mono_mass);
    list_filtered[i]->setReverseTolerance(reve_tole);
  }
  //end msThreeSetTolerance and result = list_filtered
  double ppo = peak_tole_ptr->getPpo();
  return ExtendMsPtr(new Ms<ExtendPeakPtr>(header,list_filtered, ppo));
}

} /* namespace prot */
