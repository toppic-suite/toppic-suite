#ifndef PROT_MS_HEADER_HPP_
#define PROT_MS_HEADER_HPP_

#include "activation.hpp"

namespace prot {

class MsHeader {
 public:
  MsHeader(int scan_num, int level, int charge);

  double getPrecMonoMass();

  double getPrecSpMass();

  double getPrecMonoMassMinusWater();

  std::string getScansString();

  void setScans(std::string s);

  std::string toString();

  /* get functions */
  ActivationPtr getActivationPtr() {return activation_ptr_;}

  std::string getFileName() {return file_name_;}

  int getId() {return id_;}

  int getMsLevel() {return level_;}

  std::string getTitle() {return title_;}

  double getPrecSpMz() {return prec_sp_mz_;}

  int getPrecChrg() {return prec_chrg_;}

  double getPrecMonoMz() {return prec_mono_mz_;}

  double getRetentionTime() {return retention_time_;}

  int getFirstScanNum() {return scans_[0];}

  int getPrecId() {return prec_id_;}

  /* set function */
  void setActivationPtr(ActivationPtr acti_ptr) {activation_ptr_ = acti_ptr;}

  void setFileName(std::string file_name) {file_name_ = file_name;}

  void setId(int id) {id_ = id;}

  void setTitle(std::string title) {title_ = title;}

  void setPrecSpMz(double prec_sp_mz) {prec_sp_mz_ = prec_sp_mz;}

  void setPrecChrg(int prec_chrg) {prec_chrg_ = prec_chrg;}

  void setPrecMonoMz(double prec_mono_mz) {prec_mono_mz_ = prec_mono_mz;}

  void setRetentionTime(double retention_time) {retention_time_ = retention_time;}

  void setScan(int scan_num) {scans_.push_back(scan_num);}

  void setScans(std::vector<int> scans) {scans_ = scans;}

  void setMsLevel(int level) {level_ = level;}

  void setPrecId(int prec_id) {prec_id_ = prec_id;}

 private:
  int id_ = -1;
  /** one spectrum may have several possible precursor mass */
  int prec_id_ = -1;
  /** data set name */
  std::string file_name_;

  std::string title_;
  /** a list of scans for merged spectra */
  std::vector<int> scans_;
  /** ms level */
  int level_;
  /** activation type */
  ActivationPtr activation_ptr_;
  /** retention time */
  double retention_time_;
  /** precursor m/z value in the MS1 spectrum */
  double prec_sp_mz_ = -1;
  /** computed monoisotopic precursor m/z value */
  double prec_mono_mz_ = -1;

  int prec_chrg_ = -1;
};

}

#endif
