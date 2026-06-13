// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_MS_SPEC_MS_HEADER_HPP_
#define TOPPIC_MS_SPEC_MS_HEADER_HPP_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common/base/activation.hpp"
#include "ms/spec/precursor.hpp"

namespace toppic {

class MsHeader;
using MsHeaderPtr = std::shared_ptr<MsHeader>;

class MsHeader {
 public:
  MsHeader() {}

  explicit MsHeader(XmlDOMElement element);

  // get functions for all spectra
  const std::string& getFileName() const { return file_name_; }

  int getSpecId() const { return spec_id_; }

  std::string getScansString() const;

  int getFirstScanNum() const { return scans_[0]; }

  const std::string& getTitle() const { return title_; }

  int getMsLevel() const { return level_; }

  double getRetentionTime() const { return retention_time_; }

  double getVoltage() const { return voltage_; }

  // get functions for MS/MS
  int getMsOneId() const { return ms_one_id_; }

  int getMsOneScan() const { return ms_one_scan_; }

  double getPrecTargetMz() const { return prec_target_mz_; }

  double getPrecWinBegin() const { return prec_win_begin_; }

  double getPrecWinEnd() const { return prec_win_end_; }

  ActivationPtr getActivationPtr() const { return activation_ptr_; }

  // get functions for precursor
  int getPrecNum() const { return prec_ptr_vec_.size(); }

  bool containsPrec() const { return prec_ptr_vec_.size() > 0; }

  PrecursorPtr getFirstPrecPtr() const;

  const PrecursorPtrVec& getPrecPtrVec() const { return prec_ptr_vec_; }

  int getFirstPrecId() const;

  double getFirstPrecMonoMz() const;

  int getFirstPrecCharge() const;

  double getFirstPrecInte() const;

  double getFirstPrecMonoMass() const;

  int getFirstPrecFeatureId() const;

  double getFirstPrecMonoMassMinusWater() const;

  double getFirstPrecErrorTolerance(double ppo) const;

  std::pair<int, int> getFirstPrecMonoMassMinusWaterError(double ppo,
                                                          double scale) const;

  // set functions for all spectra
  void setFileName(const std::string& file_name) { file_name_ = file_name; }

  void setSpecId(int spec_id) { spec_id_ = spec_id; }

  void setTitle(const std::string& title) { title_ = title; }

  void setScans(const std::vector<int>& scans) { scans_ = scans; }

  void setScans(const std::string& s);

  void setSingleScan(int scan_num);

  void setRetentionTime(double retention_time) {
    retention_time_ = retention_time;
  }

  void setMsLevel(int level) { level_ = level; }

  void setVoltage(double voltage) { voltage_ = voltage; }

  // set functions for MS/MS spectra
  void setMsOneId(int ms_one_id) { ms_one_id_ = ms_one_id; }

  void setMsOneScan(int ms_one_scan) { ms_one_scan_ = ms_one_scan; }

  void setPrecTargetMz(double prec_target_mz) {
    prec_target_mz_ = prec_target_mz;
  }

  void setPrecWinBegin(double prec_win_begin) {
    prec_win_begin_ = prec_win_begin;
  }

  void setPrecWinEnd(double prec_win_end) { prec_win_end_ = prec_win_end; }

  void setActivationPtr(ActivationPtr acti_ptr) {
    activation_ptr_ = std::move(acti_ptr);
  }

  // set function for precursor
  void setSinglePrecPtr(const PrecursorPtr& prec_ptr);

  void setPrecPtrVec(PrecursorPtrVec prec_ptr_vec) {
    prec_ptr_vec_ = std::move(prec_ptr_vec);
  }

  std::string toString() const;

  // Append the ms_header element under `parent` and return it.
  XmlDOMElement getHeaderXml(XmlDOMDocument* xml_doc,
                             XmlDOMElement parent) const;

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement parent) const;

  static std::string getXmlElementName() { return "ms_header"; }

  // The sorting function is used in ms2 feature assignment
  static bool cmpPrecInteDec(const MsHeaderPtr& a, const MsHeaderPtr& b);

  // Used in generating MS/MS spectra with adjusted precursor mass
  static MsHeaderPtr geneMsHeaderPtr(const MsHeaderPtr& ori_ptr,
                                     double new_prec_mass);

 private:
  // mass spec data file name
  std::string file_name_;
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
  // compensation voltage for FAIME data
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

using MsHeaderPtrVec = std::vector<MsHeaderPtr>;
using MsHeaderPtr2D = std::vector<MsHeaderPtrVec>;

}  // namespace toppic

#endif
