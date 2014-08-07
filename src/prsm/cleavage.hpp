#ifndef PROT_CLEAVAGE_HPP_
#define PROT_CLEAVAGE_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class Cleavage {
 public:
  Cleavage(int pos);

  void setPairs(PeakIonPairPtrVec pairs) {pairs_ = pairs;} 

  void setExistNIon(bool n) {exist_n_ion_ = n;};

  void setExistCIon(bool c) {exist_c_ion_ = c;};

  void setShift(double shift)  {shift_=shift;}

  void setType(const std::string &type) {type_=type;}

  void setDisplayPos(int i){display_pos_=i;}

  std::string getType(){return type_;}

  void setTrunc(const std::string &trunc) {trunc_ = trunc;}

  std::string getTrunc() {return trunc_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  int pos_;
  int display_pos_;
  bool exist_n_ion_;
  bool exist_c_ion_;
  std::string type_;
  std::string trunc_;
  double shift_;
  PeakIonPairPtrVec pairs_;
};

typedef std::shared_ptr<Cleavage> CleavagePtr;
typedef std::vector<CleavagePtr> CleavagePtrVec;

CleavagePtrVec getProteoCleavage(ProteoformPtr prot_ptr, 
                                 ExtendMsPtr ms_three_ptr,
                                 double min_mass);
} /* namespace prot */

#endif /* CLEAVAGE_HPP_ */
