#include "spec/extend_peak.hpp"

namespace prot {

ExtendPeak::ExtendPeak(DeconvPeakPtr base_peak_ptr, 
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


ExtendMsPtr createMsThreePtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, 
                             double new_prec_mass) {
  MsHeaderPtr header_ptr = getHeaderPtr(deconv_ms_ptr, new_prec_mass);

  ExtendPeakPtrVec list;
  double ext_min_mass = sp_para_ptr->getExtendMinMass();
  std::vector<double> ext_offsets = sp_para_ptr->getExtendOffsets();
  for(size_t i =0; i < deconv_ms_ptr->size(); i++){
    DeconvPeakPtr deconv_peak_ptr = deconv_ms_ptr->getPeakPtr(i);
    if(deconv_peak_ptr->getMonoMass() <= ext_min_mass) {
      double orig_mass = deconv_peak_ptr->getMonoMass();
      ExtendPeakPtr extend_peak_ptr 
          = ExtendPeakPtr(new ExtendPeak(deconv_peak_ptr,orig_mass, 1.0));
      list.push_back(extend_peak_ptr);
    }
    else{
      for(size_t j=0;j < ext_offsets.size();j++){
        double mass = deconv_peak_ptr->getMonoMass() + ext_offsets[j];
        ExtendPeakPtr extend_peak_ptr 
            = ExtendPeakPtr(new ExtendPeak(deconv_peak_ptr, mass, 1.0));
        list.push_back(extend_peak_ptr);
      }
    }
  }

  //filter extend_peak
  ExtendPeakPtrVec list_filtered;
  double min_mass = sp_para_ptr->getMinMass();
  double prec_mono_mass = header_ptr->getPrecMonoMass();
  for(size_t i =0; i < list.size();i++){
    double mass = list[i]->getPosition();
    if(mass >= min_mass && mass <= prec_mono_mass - min_mass){
      list_filtered.push_back(list[i]);
    }
  }

  // sort 
  std::sort(list_filtered.begin(),list_filtered.end(),extendPeakUp);

  //set error tolerance
  PeakTolerancePtr peak_tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  for (size_t i = 0; i < list_filtered.size();i++){
    double mass = list_filtered[i]->getBasePeakPtr()->getMonoMass();
    double ori_tole = peak_tole_ptr->compStrictErrorTole(mass);
    list_filtered[i]->setOrigTolerance(ori_tole);
    double reve_tole 
        = peak_tole_ptr->compRelaxErrorTole(mass, prec_mono_mass);
    list_filtered[i]->setReverseTolerance(reve_tole);
  }
  double ppo = peak_tole_ptr->getPpo();
  return ExtendMsPtr(new Ms<ExtendPeakPtr>(header_ptr,list_filtered, ppo));
}

inline bool massErrorUp(const std::pair<int, int> &a, const std::pair<int,int> b) {
    return a.first < b.first;
}

std::vector<std::pair<int, int>> getExtendIntMassErrorList(const ExtendMsPtrVec &ext_ms_ptr_vec, 
        double scale){
    std::vector<std::pair<int,int>> mass_errors;
    for (size_t i = 0; i < ext_ms_ptr_vec.size(); i++) {
        ExtendMsPtr ext_ms_ptr = ext_ms_ptr_vec[i];
        std::pair<int,int> last_mass_error(-1, 0);
        for(size_t j=0; j<ext_ms_ptr->size(); j++){
            int m = (int) std::round(ext_ms_ptr->getPeakPtr(j)->getPosition()*scale);
            int e = (int) std::ceil(ext_ms_ptr->getPeakPtr(j)->getOrigTolerance()*scale);
            std::pair<int,int> cur_m_e (m, e);
            if(cur_m_e.first != last_mass_error.first){
                mass_errors.push_back(cur_m_e);
                last_mass_error = cur_m_e;
            }
            else if(cur_m_e.second > last_mass_error.second){
                mass_errors.pop_back();
                mass_errors.push_back(cur_m_e);
                last_mass_error = cur_m_e;
            }
        }
    }

    std::sort(mass_errors.begin(), mass_errors.end(),massErrorUp);
    return mass_errors;
}

/*std::pair<std::vector<int>, std::vector<int>> getExtendIntMassErrorList(*/
    //ExtendMsPtr ext_ms_ptr, double scale) {
  //std::vector<int> masses;
  //std::vector<int> errors;
  //int last_mass = -1;
  //int last_error = 0;
  //for(size_t i=0;i<ext_ms_ptr->size();i++){
    //int m = (int) std::round(ext_ms_ptr->getPeakPtr(i)->getPosition()*scale);
    //int e = (int) std::ceil(ext_ms_ptr->getPeakPtr(i)->getOrigTolerance()*scale);
    //if(m != last_mass){
      //masses.push_back(m);
      //errors.push_back(e);
      //last_mass = m;
      //last_error = e;
    //}
    //else if(e>last_error){
      //errors.pop_back();
      //errors.push_back(e);
      //last_error = e;
    //}
  //}
  //std::pair<std::vector<int>, std::vector<int>> results( masses, errors);
  //return results;
/*}*/

ExtendMsPtrVec createMsThreePtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                   SpParaPtr sp_para_ptr, double new_prec_mass) {
  ExtendMsPtrVec extend_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    extend_ms_ptr_vec.push_back(
        createMsThreePtr(deconv_ms_ptr_vec[i], sp_para_ptr, new_prec_mass));
  }
  return extend_ms_ptr_vec;
}

} /* namespace prot */
