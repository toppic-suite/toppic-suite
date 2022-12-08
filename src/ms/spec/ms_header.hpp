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

namespace toppic {

class MsHeader;
typedef std::shared_ptr<MsHeader> MsHeaderPtr;

class MsHeader {
 public:
  MsHeader() {};

  MsHeader(XmlDOMElement* element);

  double getPrecMonoMass();

  double getPrecSpMass();

  double getPrecMonoMassMinusWater();

  std::pair<int,int> getPrecMonoMassMinusWaterError(double ppo, double scale);

  std::string getScansString();

  int getFirstScanNum() {return scans_[0];}

  void setScans(const std::string &s);

  std::string toString();

  // get functions 
  ActivationPtr getActivationPtr() {return activation_ptr_;}

  int getId() {return id_;}

  int getFractionId() {return fraction_id_;}

  std::string getFileName() {return file_name_;}

  int getMsLevel() {return level_;}

  std::string getTitle() {return title_;}

  int getMsOneId() {return ms_one_id_;}

  int getMsOneScan() {return ms_one_scan_;}

  double getPrecSpMz() {return prec_sp_mz_;}

  int getPrecCharge() {return prec_charge_;}

  double getPrecMonoMz();

  double getPrecTargetMz() {return prec_target_mz_;}

  double getPrecWinBegin() {return prec_target_mz_ - isolation_lower_offset_;}

  double getPrecWinEnd() {return prec_target_mz_ + isolation_upper_offset_;}

  double getRetentionTime() {return retention_time_;}

  int getPrecId() {return prec_id_;}

  double getPrecErrorTolerance(double ppo) {return getPrecMonoMass() * ppo;}

  double getPrecInte() {return prec_inte_;}

  double getTimeApex() {return time_apex_;}

  double getVoltage() {return voltage_;}

  // set function 
  void setActivationPtr(ActivationPtr acti_ptr) {activation_ptr_ = acti_ptr;}

  void setFileName(const std::string &file_name) {file_name_ = file_name;}

  void setId(int id) {id_ = id;}

  void setFractionId(int fraction_id) {fraction_id_ = fraction_id;}

  void setTitle(const std::string &title) {title_ = title;}

  void setMsOneId(int ms_one_id) {ms_one_id_ = ms_one_id;}

  void setMsOneScan(int ms_one_scan) {ms_one_scan_ = ms_one_scan;}

  void setPrecSpMz(double prec_sp_mz) {prec_sp_mz_ = prec_sp_mz;}

  void setPrecCharge(int prec_charge) {prec_charge_ = prec_charge;}

  void setPrecMonoMz(double prec_mono_mz) {prec_mono_mz_ = prec_mono_mz;}

  void setPrecTargetMz(double prec_target_mz) {prec_target_mz_ = prec_target_mz;}

  void setIsolationLowerOffset(double offset) {isolation_lower_offset_ = offset;}

  void setIsolationUpperOffset(double offset) {isolation_upper_offset_ = offset;}

  void setRetentionTime(double retention_time) {retention_time_ = retention_time;}

  void setScan(int scan_num) {scans_.push_back(scan_num);}

  void setScans(const std::vector<int> &scans) {scans_ = scans;}

  void setMsLevel(int level) {level_ = level;}

  void setPrecId(int prec_id) {prec_id_ = prec_id;}

  void setPrecInte(double inte) {prec_inte_ = inte;}

  void setVoltage(double voltage) {voltage_ = voltage;}

  XmlDOMElement* getHeaderXml(XmlDOMDocument* xml_doc);

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getXmlElementName() {return "ms_header";}

  static MsHeaderPtr geneMsHeaderPtr(MsHeaderPtr ori_ptr, double new_prec_mass);

  static bool cmpPrecInteDec(const MsHeaderPtr &a, const MsHeaderPtr &b);

 private:
  int id_ = -1;

  // a data set may have multiple fractions
  int fraction_id_ = -1;

  // mass spec data file name 
  std::string file_name_;

  // one spectrum may have several possible precursor mass */
  // precursor id 
  int prec_id_ = -1;

  std::string title_;
  // a list of scans for merged spectra 
  std::vector<int> scans_;
  // ms level 
  int level_ = 0;
  // activation type 
  ActivationPtr activation_ptr_;
  // ms1 id 
  int ms_one_id_ = -1;
  // ms1 scan number
  int ms_one_scan_ = -1;
  // retention time 
  double retention_time_ = -1;
  // precursor m/z value in the mzML file. 
  // In Thermo data, they are monoisotpic precursor m/z value  
  double prec_sp_mz_ = -1;
  // computed monoisotopic precursor m/z value 
  double prec_mono_mz_ = -1;
  // isolation window targeted m/z
  double prec_target_mz_ = -1;
  // isolation window lower offset
  double isolation_lower_offset_ = -1;
  // isolation window upper offset
  double isolation_upper_offset_ = -1;
  // precursor charge state  
  int prec_charge_ = -1;
  // precursor intensity 
  double prec_inte_ = 0;
  //rt with the highest intensity in this feature
  double time_apex_ = -1;
  //compensation voltage for FAIME data
  double voltage_ = -1;
};

typedef std::vector<MsHeaderPtr> MsHeaderPtrVec;
typedef std::vector<MsHeaderPtrVec> MsHeaderPtr2D;

}

#endif
