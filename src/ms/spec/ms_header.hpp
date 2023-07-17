//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_MS_SPEC_MS_HEADER_HPP_
#define TOPPIC_MS_SPEC_MS_HEADER_HPP_

#include "common/base/activation.hpp"
#include "ms/spec/precursor.hpp"

namespace toppic {

class MsHeader;
typedef std::shared_ptr<MsHeader> MsHeaderPtr;

class MsHeader {
 public:
  MsHeader() {};

  MsHeader(XmlDOMElement* element);

  // get functions for all spectra
  int getFractionId() {return fraction_id_;}

  std::string getFileName() {return file_name_;}

  int getSpecId() {return spec_id_;}

  std::string getScansString();

  int getFirstScanNum() {return scans_[0];}

  std::string getTitle() {return title_;}

  int getMsLevel() {return level_;}

  double getRetentionTime() {return retention_time_;}

  double getVoltage() {return voltage_;}

  // get functions for MS/MS  
  int getMsOneId() {return ms_one_id_;}

  int getMsOneScan() {return ms_one_scan_;}

  double getPrecTargetMz() {return prec_target_mz_;}

  double getPrecWinBegin() {return prec_win_begin_;}

  double getPrecWinEnd() {return prec_win_end_;}

  ActivationPtr getActivationPtr() {return activation_ptr_;}

  //get functions for precursor 
  int getPrecNum() {return prec_ptr_vec_.size();}

  bool containsPrec() {return prec_ptr_vec_.size() > 0;}

  PrecursorPtr getFirstPrecPtr();

  PrecursorPtrVec getPrecPtrVec() {return prec_ptr_vec_;}

  int getFirstPrecId();

  double getFirstPrecMonoMz();

  int getFirstPrecCharge();

  double getFirstPrecInte();

  double getFirstPrecMonoMass();

  int getFirstPrecFeatureId();

  double getFirstPrecMonoMassMinusWater();

  double getFirstPrecErrorTolerance(double ppo);

  std::pair<int,int> getFirstPrecMonoMassMinusWaterError(double ppo, double scale);

  // set functions for all spectra 
  void setFractionId(int fraction_id) {fraction_id_ = fraction_id;}

  void setFileName(const std::string &file_name) {file_name_ = file_name;}

  void setSpecId(int spec_id) {spec_id_ = spec_id;}

  void setTitle(const std::string &title) {title_ = title;}

  void setScans(const std::vector<int> &scans) {scans_ = scans;}

  void setScans(const std::string &s);

  void setSingleScan(int scan_num); 

  void setRetentionTime(double retention_time) {retention_time_ = retention_time;}

  void setMsLevel(int level) {level_ = level;}

  void setVoltage(double voltage) {voltage_ = voltage;}

  // set functions for MS/MS spectra
  void setMsOneId(int ms_one_id) {ms_one_id_ = ms_one_id;}

  void setMsOneScan(int ms_one_scan) {ms_one_scan_ = ms_one_scan;}

  void setPrecTargetMz(double prec_target_mz) {prec_target_mz_ = prec_target_mz;}

  void setPrecWinBegin(double prec_win_begin) {prec_win_begin_ = prec_win_begin;}

  void setPrecWinEnd(double prec_win_end) {prec_win_end_ = prec_win_end;}
  
  void setActivationPtr(ActivationPtr acti_ptr) {activation_ptr_ = acti_ptr;}

  // set function for precursor 
  void setSinglePrecPtr(PrecursorPtr prec_ptr); 

  void setPrecPtrVec(PrecursorPtrVec prec_ptr_vec) {prec_ptr_vec_ = prec_ptr_vec;}

  std::string toString();

  XmlDOMElement* getHeaderXml(XmlDOMDocument* xml_doc);

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getXmlElementName() {return "ms_header";}

  // The sorting function is used in ms2 feature assignment
  static bool cmpPrecInteDec(const MsHeaderPtr &a, const MsHeaderPtr &b);
  
  // Used in generating MS/MS spectra with adjusted precursor mass 
  static MsHeaderPtr geneMsHeaderPtr(MsHeaderPtr ori_ptr, double new_prec_mass);

 private:
  // mass spec data file name 
  std::string file_name_;
  // a data set may have multiple fractions
  int fraction_id_ = -1;
  // spec id
  int spec_id_ = -1;
  // mass spec title
  std::string title_;
  // a list of scans for merged spectra 
  std::vector<int> scans_;
  // ms level 
  int level_ = 0;
  // retention time 
  double retention_time_ = -1;
  //compensation voltage for FAIME data
  double voltage_ = -1;

  // information for ms/ms spectra
  // ms1 id 
  int ms_one_id_ = -1;
  // ms1 scan number
  int ms_one_scan_ = -1;
  // precursor isolation window begin
  double prec_win_begin_ = -1;
  // precusor isolation window end
  double prec_win_end_ = -1;
  // precursor isolation window targeted m/z
  double prec_target_mz_ = -1;
  // activation type 
  ActivationPtr activation_ptr_;

  // Precursor information
  PrecursorPtrVec prec_ptr_vec_;
};

typedef std::vector<MsHeaderPtr> MsHeaderPtrVec;
typedef std::vector<MsHeaderPtrVec> MsHeaderPtr2D;

}

#endif
