#ifndef PROT_CLEAVAGE_HPP_
#define PROT_CLEAVAGE_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

#define CLEAVAGE_TYPE_NORMAL "normal"
#define CLEAVAGE_TYPE_N_TRUNCATION "n_truncation"
#define CLEAVAGE_TYPE_C_TRUNCATION "c_truncation"
#define CLEAVAGE_TYPE_SEQ_START "seq_start"
#define CLEAVAGE_TYPE_SEQ_END "seq_end"

class Cleavage {
 public:
  Cleavage(int pos);

  void setPairs(PeakIonPairPtrVec pairs) {pairs_ = pairs;} 

  void setExistNIon(bool n) {exist_n_ion_ = n;};

  void setExistCIon(bool c) {exist_c_ion_ = c;};

  void setType(const std::string &type) {type_=type;}

  void setUnexpectedChange(bool u) {is_unexpected_change_ = u;}
  
  void setUnexpectedChangeColor(int color) {unexpected_change_color_ = color;}

  std::string getType(){return type_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  int pos_;
  bool exist_n_ion_;
  bool exist_c_ion_;
  std::string type_;
  PeakIonPairPtrVec pairs_;
  bool is_unexpected_change_;
  int unexpected_change_color_;
};

typedef std::shared_ptr<Cleavage> CleavagePtr;
typedef std::vector<CleavagePtr> CleavagePtrVec;

CleavagePtrVec getProteoCleavage(ProteoformPtr prot_ptr, 
                                 ExtendMsPtr ms_three_ptr,
                                 double min_mass);
} /* namespace prot */

#endif /* CLEAVAGE_HPP_ */
