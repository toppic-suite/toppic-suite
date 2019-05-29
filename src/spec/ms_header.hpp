//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_SPEC_MS_HEADER_HPP_
#define TOPPIC_SPEC_MS_HEADER_HPP_

#include "common/base/activation.hpp"

namespace toppic {

class MsHeader;
typedef std::shared_ptr<MsHeader> MsHeaderPtr;

class MsHeader {
 public:
  MsHeader() {};

  MsHeader(xercesc::DOMElement* element);

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

 // double getFeatureInte() {return feature_inte_;}

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

  double getRetentionTime() {return retention_time_;}

  int getPrecId() {return prec_id_;}

  double getErrorTolerance(double ppo) {return getPrecMonoMass() * ppo;}

  double getPrecInte() {return prec_inte_;}

  //int getFeatureId() {return feature_id_;}

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

  void setRetentionTime(double retention_time) {retention_time_ = retention_time;}

  void setScan(int scan_num) {scans_.push_back(scan_num);}

  void setScans(const std::vector<int> &scans) {scans_ = scans;}

  void setMsLevel(int level) {level_ = level;}

  void setPrecId(int prec_id) {prec_id_ = prec_id;}

  void setPrecInte(double inte) {prec_inte_ = inte;}

  //void setFeatureId(int feature_id) {feature_id_ = feature_id;}

  //void setFeatureInte(double feature_inte) {feature_inte_ = feature_inte;}

  xercesc::DOMElement* getHeaderXml(XmlDOMDocument* xml_doc);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

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
  double retention_time_ = 0.0;
  // precursor m/z value in the MS1 spectrum 
  double prec_sp_mz_ = -1;
  // computed monoisotopic precursor m/z value 
  double prec_mono_mz_ = -1;
  // precursor charge state  
  int prec_charge_ = -1;
  // precursor intensity 
  double prec_inte_ = 0;
  /*
  // feature id
  int feature_id_ = -1;
  // feature inte
  double feature_inte_ = -1;
  */
};

typedef std::vector<MsHeaderPtr> MsHeaderPtrVec;
typedef std::vector<MsHeaderPtrVec> MsHeaderPtr2D;

}

#endif
