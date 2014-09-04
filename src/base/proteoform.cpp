#include <sstream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/proteoform.hpp"

namespace prot {

Proteoform::Proteoform(DbResSeqPtr db_res_seq_ptr, 
                       ProtModPtr prot_mod_ptr,  
                       ResSeqPtr res_seq_ptr, 
                       int start_pos, int end_pos, 
                       const ChangePtrVec &change_ptr_vec) {
  db_residue_seq_ptr_ = db_res_seq_ptr;
  prot_mod_ptr_ = prot_mod_ptr;
  residue_seq_ptr_ = res_seq_ptr;
  start_pos_ = start_pos;
  end_pos_ = end_pos;
  bp_spec_ptr_ = BpSpecPtr(new BpSpec(res_seq_ptr));
  change_list_ = change_ptr_vec;
  std::sort(change_list_.begin(), change_list_.end(), compareChangeUp);
  species_id_=0;
}

Proteoform::Proteoform(xercesc::DOMElement* element, 
                       const ProteoformPtrVec &db_proteoforms){

  start_pos_ = getIntChildValue(element, "start_pos", 0);
  end_pos_ = getIntChildValue(element, "end_pos", 0);
  species_id_ = getIntChildValue(element, "species_id", 0);

  xercesc::DOMElement* db_element= getChildElement(element,"db_residue_seq",0);
  int db_seq_id = getIntChildValue(db_element, "id", 0);
  std::string db_seq_name = getChildValue(db_element, "name", 0);
  ProteoformPtr db_proteoform = db_proteoforms[db_seq_id];
  if(db_proteoform->getSeqId() != db_seq_id 
     || db_proteoform->getSeqName() != db_seq_name){
    LOG_ERROR("Sequence ID and/or name is not consistent!" << " name " 
              << db_proteoform->getSeqName() << " xml name " << db_seq_name);
    std::exit(0);
  }
  
  db_residue_seq_ptr_ = db_proteoform->getDbResSeqPtr();

  xercesc::DOMElement* mod_element= getChildElement(element,"prot_mod",0);
  std::string mod_name = getChildValue(mod_element, "name", 0);
  prot_mod_ptr_ = ProtModFactory::getBaseProtModPtrByName(mod_name);

  residue_seq_ptr_ = db_proteoform->getResSeqPtr()->getSubResidueSeq(start_pos_,end_pos_);

  if(!prot_mod_ptr_->getPtmPtr()->isEmpty()
      &&residue_seq_ptr_->getResiduePtr(0)->getPtmPtr()->isEmpty()){
    ResiduePtr mut_residue = ResidueFactory::getBaseResiduePtrByAcidPtm(residue_seq_ptr_->getResiduePtr(0)->getAcidPtr(),prot_mod_ptr_->getPtmPtr());
    ResiduePtrVec new_residue = residue_seq_ptr_->getResidues();
    new_residue[0]=mut_residue;
    residue_seq_ptr_ = ResSeqPtr(new ResidueSeq(new_residue));
  }

  bp_spec_ptr_= BpSpecPtr(new BpSpec(residue_seq_ptr_));

  xercesc::DOMElement* change_list_element= prot::getChildElement(element,"change_list",0);
  int change_len = getChildCount(change_list_element,"change");

  for(int i=0;i<change_len;i++){
    xercesc::DOMElement* change_element = getChildElement(change_list_element,"change",i);
    change_list_.push_back(ChangePtr(new Change(change_element)));
  }
}

/* get several segments without unexpected PTMs from a proteoform */
SegmentPtrVec Proteoform::getSegmentPtrVec() {
  ChangePtrVec unexpected_changes;
  double mass_shift_sum = 0;
  for (size_t i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeType() == UNEXPECTED_CHANGE) {
      unexpected_changes.push_back(change_list_[i]);
      mass_shift_sum += change_list_[i]->getMassShift();
    }
  }
  SegmentPtrVec segments;
  double n_shift = 0;
  double c_shift = mass_shift_sum;
  int left = 0;
  for (size_t i = 0; i < unexpected_changes.size(); i++) {
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
  for (size_t i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeType() == UNEXPECTED_CHANGE) {
      n++;
    }
  }
  return n;
}

ChangePtrVec Proteoform::getUnexpectedChangePtrVec(){
  ChangePtrVec un_change;
  for (size_t i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeType() == UNEXPECTED_CHANGE) {
      un_change.push_back(change_list_[i]);
    }
  }
  return un_change;
}

SemiAlignTypePtr Proteoform::getSemiAlignType() {
  int trunc_len = prot_mod_ptr_->getTruncPtr()->getTruncLen();
  //LOG_DEBUG("seq " << getProteinMatchSeq() << " trunc len " << trunc_len << " start pos " << start_pos_);
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
  double mass = getResSeqPtr()->getSeqMass();
  for(size_t i = 0;i<change_list_.size();i++){
    // only unexpected changes need to to added 
    if (change_list_[i]->getChangeType() == UNEXPECTED_CHANGE) {
      mass += change_list_[i]->getMassShift();
    }
  }
  return mass;
}

std::string Proteoform::getProteinMatchSeq(){
  std::string protein_string = db_residue_seq_ptr_->toAcidString();
  //LOG_DEBUG("protein string lenth " << protein_string.length() << " string " << protein_string);
  std::string mid_string = protein_string.substr(start_pos_, end_pos_- start_pos_ +1);
  //LOG_DEBUG("mid string lenth " << mid_string.length() << " string " << mid_string);
  std::sort(change_list_.begin(),change_list_.end(),compareChangeTypeUpPosUp);

  std::vector<std::string> break_left_strings; 
  std::vector<std::string> break_right_strings; 
  for (size_t i = 0; i < mid_string.length() + 1; i++) {
    break_left_strings.push_back("");
    break_right_strings.push_back("");
  }
  for (size_t i = 0; i < change_list_.size(); i++) {
      int left_pos = change_list_[i]->getLeftBpPos();
      break_left_strings[left_pos] = "(" + break_left_strings[left_pos];
      int right_pos = change_list_[i]->getRightBpPos();
      double shift = change_list_[i]->getMassShift();
      break_right_strings[right_pos] +=  ")";
      if (change_list_[i]->getPtmPtr() == nullptr) {
          break_right_strings[right_pos] = break_right_strings[right_pos] + "["+convertToString(shift,5)+"]";
      }
      else {
          break_right_strings[right_pos] = break_right_strings[right_pos] + "["+change_list_[i]->getPtmPtr()->getAbbrName()+"]";
      }
  }
  std::string result="";
  for (size_t i = 0; i < mid_string.length(); i++) {
    result = result + break_right_strings[i] + break_left_strings[i] + mid_string.substr(i, 1);
  }
  // last break;
  result = result + break_right_strings[mid_string.length()];
  
  std::string prefix = "";
  if(start_pos_>0){
    prefix = protein_string.substr(start_pos_-1,1);
  }
  std::string suffix = "";
  if(end_pos_< (int)protein_string.length()-1){
    suffix = protein_string.substr(end_pos_+1,1);
  }

  //LOG_DEBUG("Prefix " << prefix << " result " << result << " suffix length " << suffix.length() << " suffix " << suffix);
  return prefix+"."+result+"."+suffix;
}

void Proteoform::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("proteoform");
  std::string str = convertToString(start_pos_);
  xml_doc->addElement(element, "start_pos", str.c_str());
  str = convertToString(end_pos_);
  xml_doc->addElement(element, "end_pos", str.c_str());
  str = convertToString(species_id_);
  xml_doc->addElement(element, "species_id", str.c_str());

  db_residue_seq_ptr_->appendXml(xml_doc,element);
  prot_mod_ptr_->appendxml(xml_doc,element);
  //    residue_seq_ptr_->appendXml(xml_doc,element);
  //    bp_spec_ptr_->appendXml(xml_doc,element);
  xercesc::DOMElement* cl = xml_doc->createElement("change_list");
  for(size_t i=0;i<change_list_.size();i++){
    change_list_[i]->appendXml(xml_doc,cl);
  }
  element->appendChild(cl);
  parent->appendChild(element);
}

ProteoformPtr getDbProteoformPtr(DbResSeqPtr db_res_seq_ptr) {
  int start_pos = 0;
  int end_pos = db_res_seq_ptr->getLen() - 1;
  ChangePtrVec change_list;  
  for (int i = 0; i < db_res_seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = db_res_seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (!ptm_ptr->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, FIXED_CHANGE, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  ProtModPtr none_prot_mod_ptr = ProtModFactory::getProtModPtr_NONE();
  return ProteoformPtr(new Proteoform(db_res_seq_ptr, none_prot_mod_ptr,  
                                      db_res_seq_ptr, start_pos, end_pos, 
                                      change_list));
}

ProteoformPtr getProtModProteoform(ProteoformPtr db_form_ptr,
                                   ProtModPtr prot_mod_ptr) {
  // check if the proteoform can be truncated
  TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
  DbResSeqPtr db_res_seq_ptr = db_form_ptr->getDbResSeqPtr();  
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

ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr, 
                               int local_start, int local_end) {
  ResiduePtrVec residues;
  ResSeqPtr res_seq_ptr = proteoform_ptr->getResSeqPtr();  
  for (int i = local_start; i <= local_end; i++) {
    residues.push_back(res_seq_ptr->getResiduePtr(i));
  }
  ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
  ChangePtrVec change_list;
  ChangePtrVec ori_change_list = proteoform_ptr->getChangePtrVec();
  for (size_t i = 0; i < ori_change_list.size(); i++) {
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


ProteoformPtrVec generateProtModProteoform(const ProteoformPtrVec &ori_forms, 
                                           const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (size_t i = 0; i < ori_forms.size(); i++) {
    for (size_t j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = getProtModProteoform(ori_forms[i], prot_mods[j]);
      if (ptr.get() != nullptr) {
        new_forms.push_back(ptr);
      }
    }
  }
  return new_forms;
}

ProteoformPtrVec2D generate2DProtModProteoform(const ProteoformPtrVec &ori_forms, 
                                               const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec2D new_forms;
  for (size_t i = 0; i < ori_forms.size(); i++) {
    ProteoformPtrVec mod_forms;
    for (size_t j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = getProtModProteoform(ori_forms[i], prot_mods[j]);
      if (ptr.get() != nullptr) {
        mod_forms.push_back(ptr);
      }
    }
    new_forms.push_back(mod_forms);
  }
  return new_forms;
}

ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts;
  ResiduePtrVec residue_list;
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
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
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(), 
                                            residue_list[i]->getPtmPtr(),
                                            counts[i]/sum));
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}


ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list, 
                              const ProteoformPtrVec &prot_mod_forms) {
  std::vector<double> counts(residue_list.size(), 0.0);
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
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
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(), 
                                            residue_list[i]->getPtmPtr(),
                                            counts[i]/sum));
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

void Proteoform::addUnexpectedChangePtrVec(const ChangePtrVec &changes) {
  for (size_t i = 0; i < changes.size(); i++) {
    if (changes[i]->getChangeType() == UNEXPECTED_CHANGE) {
      change_list_.push_back(changes[i]);
    }
  }
}

bool isSamePeptideAndMass(ProteoformPtr a, ProteoformPtr b, double ppo) {  
  if(a->getDbResSeqPtr()->getId() != b->getDbResSeqPtr()->getId()){
    return false;
  }
  if(a->getStartPos() != b->getStartPos()){
    return false;
  }
  if(a->getEndPos() != b->getEndPos()){
    return false;
  }
  double thresh = a->getResSeqPtr()->getSeqMass() * ppo;
  if(std::abs(a->getResSeqPtr()->getSeqMass()
              -b->getResSeqPtr()->getSeqMass())> thresh){
    return false;
  }
  return true;
}

bool isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b,
                                   double ppo){
  if(!isSamePeptideAndMass(a,b,ppo)){
    return false;
  }
  if(a->getChangePtrVec().size() != b->getChangePtrVec().size()){
    return false;
  }
  double shift_tolerance = a->getResSeqPtr()->getSeqMass()*ppo;
  // sort changes 
  ChangePtrVec a_change_vec = a->getChangePtrVec();
  ChangePtrVec b_change_vec = b->getChangePtrVec();
  std::sort(a_change_vec.begin(),a_change_vec.end(),compareChangeUp);
  std::sort(b_change_vec.begin(),b_change_vec.end(),compareChangeUp);
  for(size_t i=0;i< a->getChangePtrVec().size();i++){
    ChangePtr ac = a_change_vec[i];
    ChangePtr bc = b_change_vec[i];
    if(ac->getRightBpPos() <= bc->getLeftBpPos() || bc->getRightBpPos() <= ac->getLeftBpPos()){
      return false;
    }
    if(abs(ac->getMassShift()-bc->getMassShift() > shift_tolerance)){
      return false;
    }
  }
  return true;
}

ProteoformPtrVec2D getProteoBlocks(const ProteoformPtrVec &proteo_ptrs, int db_block_size) {
  size_t start_idx = 0;
  size_t proteo_idx =0;
  int block_len=0;
  ProteoformPtrVec2D proteo_blocks;
  while(proteo_idx < proteo_ptrs.size()){
    int proteo_len = proteo_ptrs[proteo_idx]->getResSeqPtr()->getLen();
    if(block_len + proteo_len < db_block_size) {
      block_len = block_len + proteo_len;
      proteo_idx++;
    }
    else{
      size_t end_idx = proteo_idx;
      ProteoformPtrVec proteo_in_block;
      for(size_t i = start_idx; i<=end_idx; i++){
        proteo_in_block.push_back(proteo_ptrs[i]);
      }
      proteo_blocks.push_back(proteo_in_block);
      start_idx = end_idx +1;
      proteo_idx = end_idx +1;
      block_len = 0;
    }
  }
  /* last block */
  if (start_idx < proteo_ptrs.size()) {
    ProteoformPtrVec proteo_in_block;
    for(size_t i = start_idx; i < proteo_ptrs.size(); i++){
      proteo_in_block.push_back(proteo_ptrs[i]);
    }
    proteo_blocks.push_back(proteo_in_block);
  }
  return proteo_blocks;
}

} /* namespace prot */

