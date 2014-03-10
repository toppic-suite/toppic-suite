#include <sstream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/proteoform.hpp"

namespace prot {

Proteoform::Proteoform(DbResSeqPtr db_res_seq_ptr, ProtModPtr prot_mod_ptr,  
                       ResSeqPtr res_seq_ptr, int start_pos, int end_pos, 
                       ChangePtrVec change_list) {
  db_residue_seq_ptr_ = db_res_seq_ptr;
  prot_mod_ptr_ = prot_mod_ptr;
  residue_seq_ptr_ = res_seq_ptr;
  start_pos_ = start_pos;
  end_pos_ = end_pos;
  LOG_TRACE( "start bp spec generation");
  bp_spec_ptr_ = BpSpecPtr(new BpSpec(res_seq_ptr));
  change_list_ = change_list;
  std::sort(change_list.begin(), change_list.end(), compareChangeUp);
  species_id_=0;
}

Proteoform::Proteoform(xercesc::DOMElement* element,ProteoformPtrVec proteoforms){

  start_pos_ = getIntChildValue(element, "start_pos", 0);
  end_pos_ = getIntChildValue(element, "end_pos", 0);

  xercesc::DOMElement* db_element= prot::getChildElement(element,"db_residue_seq",0);
  int db_seq_id = getIntChildValue(db_element, "id", 0);
  std::string db_seq_name = getChildValue(db_element, "name", 0);
  ProteoformPtr seq = proteoforms[db_seq_id];
  if(seq->getSeqId() != db_seq_id || seq->getName().compare(db_seq_name)!=0){
    std::cout<< "Sequence ID and/or name is not consistent!" << std::endl;
    std::exit(0);
  }
  db_residue_seq_ptr_ = seq->getDbResSeqPtr();

  xercesc::DOMElement* mod_element= prot::getChildElement(element,"prot_mod",0);
  std::string mod_name = getChildValue(mod_element, "name", 0);
  //    double mod_prot_shift = getDoubleChildValue(mod_element, "prot_shift", 0);
  //    double mod_pep_shift = getDoubleChildValue(mod_element, "pep_shift", 0);
  //    prot_mod_ptr_ = getProtModPtrByName(basedata->getProtModPtrVec(),mod_name);
  prot_mod_ptr_ = ProtModFactory::getBaseProtModPtrByName(mod_name);

//  xercesc::DOMElement* res_seq_element= getChildElement(element,"residue_seq",0);
//  int res_len = getChildCount(res_seq_element,"residue");
//  ResiduePtrVec residues;
//  for(int i=0;i<res_len;i++){
//    xercesc::DOMElement* res_element= getChildElement(res_seq_element,"residue",i);
//    std::string acid_name
//        = getChildValue(getChildElement(res_element,"amino_acid",0),"name",0);
//    std::string ptm_name
//        = getChildValue(getChildElement(res_element,"modification",0),"abbr_name",0);
//    residues.push_back(ResiduePtr(new Residue(acid_name,ptm_name)));
//  }
//  residue_seq_ptr_ = ResSeqPtr(new ResidueSeq(residues));
  residue_seq_ptr_ = proteoforms[db_residue_seq_ptr_->getId()]->getResSeqPtr();
  bp_spec_ptr_= BpSpecPtr(new BpSpec(residue_seq_ptr_));

  xercesc::DOMElement* change_list_element= prot::getChildElement(element,"change_list",0);
  int change_len = getChildCount(change_list_element,"change");

  for(int i=0;i<change_len;i++){
    xercesc::DOMElement* change_element
        = prot::getChildElement(change_list_element,"change",i);
    int left_bp_pos = getIntChildValue(change_element,"left_bp_pos",0);
    int right_bp_pos = getIntChildValue(change_element,"right_bp_pos",0);
    int change_type = getIntChildValue(change_element,"change_type",0);
    double mass_shift = getDoubleChildValue(change_element,"mass_shift",0);
    int ptm_count = getChildCount(change_element,"modification");
    PtmPtr change_ptm = nullptr;
    if(ptm_count!=0){
      xercesc::DOMElement* ptm_element
          = getChildElement(change_element,"modification",i);
      change_ptm 
          = PtmFactory::getBasePtmPtrByAbbrName(getChildValue(ptm_element,"abbr_name",0));
    }
    change_list_.push_back(
        ChangePtr(new Change(left_bp_pos,right_bp_pos,change_type,mass_shift,change_ptm)));
  }
  species_id_=0;
}

SegmentPtrVec Proteoform::getSegmentPtrVec() {
  ChangePtrVec unexpected_changes;
  double mass_shift_sum = 0;
  for (unsigned int i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeType() == UNEXPECTED_CHANGE) {
      unexpected_changes.push_back(change_list_[i]);
      mass_shift_sum += change_list_[i]->getMassShift();
    }
  }
  SegmentPtrVec segments;
  double n_shift = 0;
  double c_shift = mass_shift_sum;
  int left = 0;
  for (unsigned int i = 0; i < unexpected_changes.size(); i++) {
    int right = unexpected_changes[i]->getLeftBpPos();
    SegmentPtr segment_ptr = SegmentPtr(
        new Segment(left, right, n_shift, c_shift)); 
    segments.push_back(segment_ptr);
    left = unexpected_changes[i]->getRightBpPos();
    n_shift = n_shift + unexpected_changes[i]->getMassShift();
    c_shift = c_shift - unexpected_changes[i]->getMassShift();
  }
  int right = residue_seq_ptr_->getLen();
  SegmentPtr segment_ptr = SegmentPtr(
      new Segment(left, right, n_shift, c_shift)); 
  segments.push_back(segment_ptr);
  return segments;
}

std::string Proteoform::toString() {
  std::stringstream s;
  s<< "Begin pos: " << start_pos_ << std::endl;
  s<< "End pos: " << end_pos_ << std::endl;
  s<< "String: " << residue_seq_ptr_->toString();
  return s.str();
}

int Proteoform::getUnexpectedChangeNum() {
  int n = 0;
  for (unsigned int i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeType() == UNEXPECTED_CHANGE) {
      n++;
    }
  }
  return n;
}

SemiAlignTypePtr Proteoform::getSemiAlignType() {
  int trunc_len = prot_mod_ptr_->getTruncPtr()->getTruncLen();
  bool is_prefix = false;
  if (start_pos_ == trunc_len) {
    is_prefix = true;
  }

  bool is_suffix = false;
  if (end_pos_ == db_residue_seq_ptr_->getLen() - 1) {
    is_suffix = true;
  }

  if (is_prefix) {
    if (is_suffix) {
      return SemiAlignTypeFactory::getCompletePtr();
    } else {
      return SemiAlignTypeFactory::getPrefixPtr();
    }
  } else {
    if (is_suffix) {
      return SemiAlignTypeFactory::getSuffixPtr();
    } else {
      return SemiAlignTypeFactory::getInternalPtr();
    }
  }
}

double Proteoform::getMass(){
  std::vector<double> ext_b_mass 
      = bp_spec_ptr_->getBreakPointMasses(IonTypeFactory::getIonTypePtr_B());
  double mass = ext_b_mass[end_pos_]-ext_b_mass[start_pos_]+ MassConstant::getWaterMass();
  for(unsigned int i = 0;i<change_list_.size();i++){
    mass += change_list_[i]->getMassShift();
  }
  return mass;
}

bool Proteoform::isAcetylation(){
  if(prot_mod_ptr_->getPtmPtr()->getAbbrName().compare("Acetylation") == 0){
    return true;
  }
  else{
    return false;
  }
}

std::string Proteoform::getProteinMatchSeq(){
  std::string result="";
  std::string protein_string = residue_seq_ptr_->toString();
  std::string mid_string = protein_string.substr(start_pos_,end_pos_);
  int mid_start=0;
  std::sort(change_list_.begin(),change_list_.end(),compareChangeUp);
  for(unsigned int i=0;i<change_list_.size();i++){
    result += mid_string.substr(mid_start,change_list_[i]->getLeftBpPos()-mid_start);
    mid_start=change_list_[i]->getLeftBpPos();
    result += "(";
    result += mid_string.substr(mid_start,change_list_[i]->getRightBpPos()-mid_start);
    result += ")";
    result += "["+convertToString(change_list_[i]->getMassShift(),5)+"]";
    mid_start = change_list_[i]->getRightBpPos();
  }
  result += mid_string.substr(mid_start,end_pos_-start_pos_-mid_start);

  std::string prefix = "";
  if(start_pos_>0){
    prefix = protein_string.substr(start_pos_-1,1);
  }
  std::string suffix = "";

  if(end_pos_< (int)protein_string.length()-2){
    suffix = protein_string.substr(end_pos_,1);
  }

  return prefix+"."+result+"."+suffix;
}

void Proteoform::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
    xercesc::DOMElement* element = xml_doc->createElement("proteoform");
    std::string str = convertToString(start_pos_);
    xml_doc->addElement(element, "start_pos", str.c_str());
    str = convertToString(end_pos_);
    xml_doc->addElement(element, "end_pos", str.c_str());
    db_residue_seq_ptr_->appendXml(xml_doc,element);
    prot_mod_ptr_->appendxml(xml_doc,element);
//    residue_seq_ptr_->appendXml(xml_doc,element);
//    bp_spec_ptr_->appendXml(xml_doc,element);
    xercesc::DOMElement* cl = xml_doc->createElement("change_list");
    for(unsigned int i=0;i<change_list_.size();i++){
        change_list_[i]->appendXml(xml_doc,cl);
    }
    element->appendChild(cl);
    parent->appendChild(element);
}

ProteoformPtr getDbProteoformPtr(DbResSeqPtr db_res_seq_ptr, 
                                 ProtModPtr prot_mod_ptr) {
  int start_pos = 0;
  int end_pos = db_res_seq_ptr->getLen() - 1;
  //LOG_DEBUG("raw protein sequence name " << db_res_seq_ptr->getName() 
  //<< " len " << db_res_seq_ptr->getLen());
  ChangePtrVec change_list;  
  for (int i = 0; i < db_res_seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = db_res_seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (!ptm_ptr->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, FIXED_CHANGE, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  return ProteoformPtr(new Proteoform(db_res_seq_ptr, prot_mod_ptr,  
                                      db_res_seq_ptr, start_pos, end_pos, 
                                      change_list));
}

ProteoformPtr getProtModProteoform(ProteoformPtr raw_form_ptr,
                                   ProtModPtr prot_mod_ptr) {
  // check if the proteoform can be truncated
  //LOG_DEBUG("Prot mod " << prot_mod_ptr->getName());
  TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
  DbResSeqPtr db_res_seq_ptr = raw_form_ptr->getDbResSeqPtr();  
  bool valid_trunc = trunc_ptr->isValidTrunc(db_res_seq_ptr);
  if (!valid_trunc) {
    //LOG_DEBUG("NO valid trunc");
    return ProteoformPtr(nullptr);
  }
  // first residue might be acetylated 
  int start = trunc_ptr->getTruncLen();
  ResiduePtrVec residues;
  ResiduePtr residue = db_res_seq_ptr->getResiduePtr(start);
  PtmPtr ori_ptm = residue->getPtmPtr();
  PtmPtr ptm = prot_mod_ptr->getPtmPtr();
  ChangePtrVec change_list;
  if (ptm->isEmpty() || !ori_ptm->isEmpty()) {
    residues.push_back(residue);
    if (!ori_ptm->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(0, 1, FIXED_CHANGE, ori_ptm->getMonoMass(), ori_ptm));
      change_list.push_back(change_ptr);
    }
  }
  else {
    AcidPtr acid = residue->getAcidPtr();
    ResiduePtr mut_residue = ResidueFactory::getBaseResiduePtrByAcidPtm(acid, ptm);
    if (mut_residue.get() == nullptr) {
      LOG_ERROR( "Proteoform:: residue not found");
      throw("Residue not found");
    }
    residues.push_back(mut_residue);
    change_list.push_back(ChangePtr(new Change(0,1, PROTEIN_VARIABLE_CHANGE, 
                                               ptm->getMonoMass(), ptm)));
  }
  // add all other residues
  for (int i = start + 1; i < db_res_seq_ptr->getLen(); i++) {
    residues.push_back(db_res_seq_ptr->getResiduePtr(i));
  }
  ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
  // start from 1 since ptm on the first residue has been added to the change
  // list
  for (int i = 1; i < seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (!ptm_ptr->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, FIXED_CHANGE, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  //LOG_DEBUG("mod protein sequence name " << db_res_seq_ptr->getName() 
  //<< " len " << db_res_seq_ptr->getLen());
  return ProteoformPtr(
      new Proteoform(db_res_seq_ptr, prot_mod_ptr, seq_ptr, start, 
                     db_res_seq_ptr->getLen()-1, change_list));
}

ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr, int local_start, 
                               int local_end) {
  ResiduePtrVec residues;
  ResSeqPtr res_seq_ptr = proteoform_ptr->getResSeqPtr();  
  for (int i = local_start; i <= local_end; i++) {
    residues.push_back(res_seq_ptr->getResiduePtr(i));
  }
  ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
  ChangePtrVec change_list;
  ChangePtrVec ori_change_list = proteoform_ptr->getChangePtrVec();
  for (unsigned int i = 0; i < ori_change_list.size(); i++) {
    if (ori_change_list[i]->getLeftBpPos() >= local_start 
        && ori_change_list[i]->getRightBpPos() <= local_end + 1) {
      ChangePtr change_ptr = ChangePtr(new Change(*ori_change_list[i], local_start));
      change_list.push_back(change_ptr);
    }
  }
  DbResSeqPtr db_res_seq_ptr = proteoform_ptr->getDbResSeqPtr();
  ProtModPtr prot_mod_ptr = proteoform_ptr->getProtModPtr();
  return ProteoformPtr(
      new Proteoform(db_res_seq_ptr, prot_mod_ptr, seq_ptr, 
                     local_start + proteoform_ptr->getStartPos(), 
                     local_end + proteoform_ptr->getStartPos(), change_list));
}


ProteoformPtrVec generateProtModProteoform(ProteoformPtrVec &ori_forms, 
                                           ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (unsigned int i = 0; i < ori_forms.size(); i++) {
    for (unsigned int j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = getProtModProteoform(ori_forms[i], prot_mods[j]);
      if (ptr.get() != nullptr) {
        new_forms.push_back(ptr);
      }
    }
  }
  return new_forms;
}

ResFreqPtrVec compNTermResidueFreq(ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts;
  ResiduePtrVec residue_list;
  for (unsigned int i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();    
    if (seq_ptr->getLen() >= 1) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
      int pos = findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found 
        counts[pos] = counts[pos]+1;
      }
      else {
        residue_list.push_back(res_ptr);
        counts.push_back(1);
      }
    }
  }

  double sum = 0;
  for (unsigned int i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (unsigned int i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(), 
                                            residue_list[i]->getPtmPtr(),
                                            counts[i]/sum));
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}


ResFreqPtrVec compResidueFreq(ResiduePtrVec &residue_list, 
                              ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts(residue_list.size(), 0.0);
  for (unsigned int i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();    
    for (int j = 0; j < seq_ptr->getLen(); j++) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(j);
      int pos = findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found 
        counts[pos] = counts[pos]+1;
      }
    }
  }

  double sum = 0;
  for (unsigned int i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (unsigned int i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(), 
                                            residue_list[i]->getPtmPtr(),
                                            counts[i]/sum));
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

void Proteoform::addUnexpectedChangePtrVec(ChangePtrVec &changes) {
  for (unsigned int i = 0; i < changes.size(); i++) {
    change_list_.push_back(changes[i]);
  }
}

bool isSamePeptideAndMass(ProteoformPtr proteoform,ProteoformPtr another_proteoform,double ppo){
  double thresh = proteoform->getBpSpecPtr()->getResSeqMass()*ppo;
  if(proteoform->getDbResSeqPtr()->getId() != another_proteoform->getDbResSeqPtr()->getId()){
    return false;
  }
  if(proteoform->getStartPos() != another_proteoform->getStartPos()){
    return false;
  }
  if(proteoform->getEndPos() != another_proteoform->getEndPos()){
    return false;
  }
  if(std::abs(proteoform->getBpSpecPtr()->getResSeqMass()
              -another_proteoform->getBpSpecPtr()->getResSeqMass())> thresh){
    return false;
  }
  return true;
}

bool isStrictCompatiablePtmSpecies(ProteoformPtr a,ProteoformPtr b,double ppo){
  if(!isSamePeptideAndMass(a,b,ppo)){
    return false;
  }
  if(a->getChangePtrVec().size() != b->getChangePtrVec().size()){
    return false;
  }
  double shift_tolerance = a->getBpSpecPtr()->getResSeqMass()*ppo;
  for(unsigned int i=0;i< a->getChangePtrVec().size();i++){
    ChangePtr ac = a->getChangePtrVec()[i];
    ChangePtr bc = b->getChangePtrVec()[i];
    if(ac->getRightBpPos() <= bc->getLeftBpPos() || bc->getRightBpPos() <= ac->getLeftBpPos()){
      return false;
    }
    if(abs(ac->getMassShift()-bc->getMassShift() > shift_tolerance)){
      return false;
    }
  }
  return true;
}

} /* namespace prot */

