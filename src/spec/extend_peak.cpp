/*
 * extend_peak.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include "spec/extend_peak.hpp"

namespace prot {

ExtendPeak::ExtendPeak(DeconvPeakPtr base_peak,double mono_mass,double score):Peak(mono_mass,1.0){
	base_peak_ = base_peak;
	mono_mass_ = mono_mass;
	score_ = score;
	orig_tolerance_ = 0.0;
	reverse_tolerance_ = 0.0;
}

void ExtendPeak::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
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
	base_peak_->appendXml(xml_doc,element);
	parent->appendChild(element);
}

//must deleted the ms after finished using
//ExtendMsPtr getMsThree(DeconvMsPtr deconv_ms,double delta,SpParaPtr sp_para);
ExtendMsPtr getMsThree(DeconvMsPtr deconv_ms,double delta,SpParaPtr sp_para){
	MsHeaderPtr header = getDeltaHeaderPtr(deconv_ms,delta);

	//private function getSpThreeExtendPeak in factory
	ExtendPeakPtrVec list;
	for(unsigned int i =0;i<deconv_ms->size();i++){
		if(deconv_ms->getPeakPtr(i)->getMonoMass() <= sp_para->getExtendSpPara()->extend_min_mass_){
			double orig_mass = deconv_ms->getPeakPtr(i)->getMonoMass();
			ExtendPeakPtr new_peak = ExtendPeakPtr(new ExtendPeak(deconv_ms->getPeakPtr(i),orig_mass,1));
			list.push_back(new_peak);
		}
		else{
			for(unsigned int j=0;j<sp_para->getExtendSpPara()->ext_offsets_.size();j++){
				double mass = deconv_ms->getPeakPtr(i)->getMonoMass() +
						sp_para->getExtendSpPara()->ext_offsets_[j];
				ExtendPeakPtr new_peak = ExtendPeakPtr(new ExtendPeak(deconv_ms->getPeakPtr(i),mass,1));
				list.push_back(new_peak);
			}
		}
	}
	//filter extend_peak
	ExtendPeakPtrVec list_filtered;
	double pre_mono_mass = header->getPrecMonoMass();
	for(unsigned int i =0; i < list.size();i++){
		double mass = list[i]->getPosition();
		if(mass >= sp_para->getMinMass() && mass <= pre_mono_mass - sp_para->getMinMass()){
			list_filtered.push_back(list[i]);
		}
	}
	//end filter extend_peak function and result = list_filtered

	std::sort(list_filtered.begin(),list_filtered.end(),extendpeak_up);
	//end getSpThreeExtendPeak and result = list_filtered;

	//function msThreeSetTolerance
	for (unsigned int i = 0; i < list_filtered.size();i++){
		list_filtered[i]->setOrigTolerance(
				sp_para->getPeakTolerance()->
        compStrictErrorTole(list_filtered[i]->getBasePeakPtr()->getMonoMass()));
		list_filtered[i]->setReverseTolerance(
				sp_para->getPeakTolerance()->compRelaxErrorTole(
						list_filtered[i]->getBasePeakPtr()->getMonoMass(),pre_mono_mass));
	}
	//end msThreeSetTolerance and result = list_filtered
  double ppo = sp_para->getPeakTolerance()->getPpo();
	 return ExtendMsPtr(new Ms<ExtendPeakPtr>(header,list_filtered, ppo));
}

} /* namespace prot */
