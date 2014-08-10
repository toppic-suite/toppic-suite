#ifndef PROT_PRSM_HPP_
#define PROT_PRSM_HPP_

#include <string>

#include "base/extreme_value.hpp"
#include "base/proteoform.hpp"
#include "base/anno_residue.hpp"
#include "spec/deconv_peak.hpp"
#include "spec/extend_peak.hpp"
#include "spec/sp_para.hpp"
#include "prsm/cleavage.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class Prsm {
 public:
  Prsm(ProteoformPtr proteoform_ptr, DeconvMsPtr deconv_ms_ptr, 
       double adjusted_prec_mass, double calibration, SpParaPtr sp_para_ptr);

  Prsm(xercesc::DOMElement* element,ProteoformPtrVec proteoforms);

  double getAdjustedPrecMass() {return adjusted_prec_mass_;}

  ProteoformPtr getProteoformPtr() {return proteoform_ptr_;}

  void setProteoformPtr(ProteoformPtr proteoform) {proteoform_ptr_=proteoform;}

  double getCalibration() {return calibration_;}

  DeconvMsPtr getDeconvMsPtr() {return deconv_ms_ptr_;}

  void setDeconvMsPtr(DeconvMsPtr ms) {deconv_ms_ptr_=ms;}

  double getMatchPeakNum() {return match_peak_num_;}

  double getEValue(); 

  double getPValue();

  double getOneProtProb();

  double getFdr() {return fdr_;}

  int getId() {return prsm_id_;}

  double getMatchFragNum() {return match_fragment_num_;}

  int getSpectrumId() {return spectrum_id_;}

  std::string getSpectrumScan() {return spectrum_scan_;}

  double getOriPrecMass() {return ori_prec_mass_;}

  int getPrecursorId() {return precursor_id_;}

  ExtendMsPtr getRefineMs() {return refine_ms_three_;}

  bool isMatchMs(MsHeaderPtr header_ptr);

  void setRefineMs(ExtendMsPtr refine_ms_three){refine_ms_three_ = refine_ms_three;}

  ExtremeValuePtr getProbPtr() {return prob_ptr_;} 

  void setProbPtr(ExtremeValuePtr prob_ptr) {prob_ptr_ = prob_ptr;}

  void setFdr(double fdr) {fdr_ = fdr;}

  void setId(int id) {prsm_id_ = id;}

  void setSpectrumId(int spectrum_id) {spectrum_id_ = spectrum_id;}

  void setSpectrumScan(std::string spectrum_scan) {spectrum_scan_ = spectrum_scan;}

  void setOriPrecMass(double prec_mass) {ori_prec_mass_ = prec_mass;}

  void setPrecurorId(int precursor_id) {precursor_id_ = precursor_id;}

  //int getMinMass(){return sp_para_ptr_->getMinMass();}

  xercesc::DOMElement* toXmlElement(XmlDOMDocument* xml_doc);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  int prsm_id_ = -1;
  /* spectrum information */
  int spectrum_id_;

  std::string spectrum_scan_;

  int precursor_id_;

  double ori_prec_mass_;
  /* adjusted precursor mass */
  double adjusted_prec_mass_;
  /* calibration */
  double calibration_;

  /* protein sequence */
  ProteoformPtr proteoform_ptr_;

  ExtremeValuePtr prob_ptr_;

  double fdr_ = -1;

  /* The following are not saved in xml */
  DeconvMsPtr deconv_ms_ptr_;
  /* adjusted extended msThree, used for matching ions and peaks */
  ExtendMsPtr refine_ms_three_;

  double match_peak_num_ = 0;
  double match_fragment_num_ = 0;

  //SpParaPtr sp_para_ptr_;

  void init(SpParaPtr sp_para_ptr);
  void initScores(SpParaPtr sp_para_ptr);
};

typedef std::shared_ptr<Prsm> PrsmPtr;
typedef std::vector<PrsmPtr> PrsmPtrVec;
typedef std::vector<PrsmPtrVec> PrsmPtrVec2D;
typedef std::vector<PrsmPtrVec2D> PrsmPtrVec3D;

inline bool prsmEValueUp(const PrsmPtr &a, const PrsmPtr &b) {
  return a->getEValue() < b->getEValue();
}

inline bool prsmEValueDown(const PrsmPtr &a, const PrsmPtr &b) {
  return a->getEValue() > b->getEValue();
}

inline bool prsmMatchFragmentDown(const PrsmPtr &a, const PrsmPtr &b) {
  return a->getMatchFragNum() > b->getMatchFragNum();
}

// sort by the order of protein sequences in the database 
inline bool prsmProteoformIdUp(PrsmPtr p1, PrsmPtr p2) {
  return p1->getProteoformPtr()->getDbResSeqPtr()->getId()
      < p2->getProteoformPtr()->getDbResSeqPtr()->getId();
}

// sort by the order of spectrum id, the precursor id
inline bool prsmSpectrumIdUpPrecursorIdUp(const PrsmPtr &a, const PrsmPtr &b){
  if(a->getSpectrumId() < b->getSpectrumId()){
    return true;
  }
  else if(a->getSpectrumId() > b->getSpectrumId()){
    return false;
  }
  else{
    if(a->getPrecursorId() < b->getPrecursorId()){
      return true;
    }
    return false;
  }
}

// sort by number of matched fragment ions, then start position 
inline bool prsmCompPreMatchFragDown(const PrsmPtr &a, const PrsmPtr &b) {
  if(a->getMatchFragNum() > b->getMatchFragNum()){
    return true;
  }
  else if(a->getMatchFragNum() == b->getMatchFragNum()){
    return a->getProteoformPtr()->getStartPos()<b->getProteoformPtr()->getStartPos();
  }
  return false;
}

// sort by spectrum id then match ions
inline bool prsmSpectrumIdUpMatchFragDown(const PrsmPtr &a, const PrsmPtr &b){
  if(a->getSpectrumId() < b->getSpectrumId()){
    return true;
  }
  else if(a->getSpectrumId() > b->getSpectrumId()){
    return false;
  }
  else{
    if(a->getMatchFragNum() > b->getMatchFragNum()){
      return true;
    }
    return false;
  }
}

// sort by spectrum id then evalue
inline bool prsmSpectrumIdUpEvalueUp(const PrsmPtr &a, const PrsmPtr &b){
  if(a->getSpectrumId() < b->getSpectrumId()){
    return true;
  }
  else if(a->getSpectrumId() > b->getSpectrumId()){
    return false;
  }
  else{
    if(a->getEValue() < b->getEValue()){
      return true;
    }
    return false;
  }
}

PrsmPtrVec readPrsm(const std::string &file_name, const ProteoformPtrVec &proteoforms);

void filterPrsms(const PrsmPtrVec &prsms, MsHeaderPtr header_ptr, PrsmPtrVec &sele_prsms); 

void addSpectrumPtrsToPrsms(PrsmPtrVec &prsms, PrsmParaPtr prsm_para_ptr);

}
#endif

