/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_RESIDUE_HPP_
#define PROT_RESIDUE_HPP_

#include <string>
#include <memory>

#include "base/acid.hpp"
#include "base/ptm.hpp"

namespace prot {

class Residue {
 public:
  Residue(const AcidPtr &acid_ptr, const PtmPtr &ptm_ptr); 

  Residue(const std::string &acid_name, const std::string &abbr_name);

  /** Get amino acid. */
  AcidPtr getAcidPtr() {return acid_ptr_; }

  /** Get residue mass. */
  double getMass() { return mass_; }

  /** Get post-translational modification. */
  PtmPtr getPtmPtr() { return ptm_ptr_; }

  void setPos(int pos) { pos_ = pos; }

  void setDisplayPos(int pos) {display_pos_ = pos;}

  void setType(std::string type) {type_ = type;}

  void setType(std::string type) {type_=type;}

  void setShiftStyle(std::string shift_type) {shift_style_=shift_type;}

  std::string getType(){return type_;}

  void setIsModifyed(bool isMod) {is_modified_ = isMod;}

  bool isModified() {return is_modified_;}

  void setShift(double shift) {shift_ = shift;}

  double getShift(){return shift_;}

  void setExpected(bool e) { is_expected_ = e;}

  /**
   * Checks if the residue contains the same amino acid and ptm.
   */
  bool isSame(const AcidPtr &acid_ptr, 
              const PtmPtr &ptm_ptr) {
    return acid_ptr_.get() == acid_ptr.get() && ptm_ptr_.get() == ptm_ptr.get();
  }

  /** Get string representation */
  std::string toString(const std::string &delim_bgn, 
                       const std::string &delim_end);

  std::string toString() {return toString("[", "]");}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  void appendViewXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  /** amino acid */
  AcidPtr acid_ptr_;
  /** post-translational modification */
  PtmPtr ptm_ptr_;
  /** residue mass */
  double mass_;

  //residues pos
  int pos_=0;
  //residues type
  std::string type_="normal";
  //display pos shown in match seq
  int display_pos_=0;
  //is modified
  bool is_modified_ = false;
  //is expected
  bool is_expected_ = false;
  //expected shift type
  std::string shift_style_;
  //modify mass
  double shift_=0;
};

typedef std::shared_ptr<Residue> ResiduePtr;
typedef std::vector<ResiduePtr> ResiduePtrVec;

ResiduePtr getResiduePtrByAcid(const ResiduePtrVec &residue_list,
                               const AcidPtr &acid_ptr);

int findResidue(const ResiduePtrVec &residue_list, const ResiduePtr &residue_ptr);

ResiduePtrVec convertAcidToResidueSeq(const ResiduePtrVec &residue_list,
                                      const AcidPtrVec &acid_list);
/* residue factory */
class ResidueFactory {
 public:
  static void initFactory(const std::string &file_name);

  static ResiduePtrVec& getBaseResiduePtrVec() {return residue_ptr_vec_;}
  
  static ResiduePtr getBaseResiduePtrByAcidPtm(const AcidPtr &acid_ptr, 
                                               const PtmPtr &ptm_ptr);
  
  static ResiduePtr addBaseResidue(const AcidPtr &acid_ptr, 
                                   const PtmPtr &ptm_ptr);
  
  static ResiduePtrVec getResiduePtrVecInstance(const std::string &file_name);

 private:
  static ResiduePtrVec residue_ptr_vec_;
};

}
#endif
