//
// Created by abbash on 8/22/22.
//

#ifndef TOPPIC_MS_MS_MAP_MS_MAP_ROW_HEADER_HPP
#define TOPPIC_MS_MS_MAP_MS_MAP_ROW_HEADER_HPP

#include <memory>
#include <vector>

namespace toppic {

class MsMapRowHeader {
 public:
  MsMapRowHeader(int row_id, int spec_id, int scan_num, double rt);

  int getRowId() {return row_id_;}

  int getRawSpecId() const { return raw_spec_id_; }
  void setRawSpecId(int raw_spec_id) { raw_spec_id_ = raw_spec_id; }

  int getScanNum() const { return scan_num_; }
  void setScanNum(int scan_num) { scan_num_ = scan_num; }

  double getRt() const { return rt_; }
  void setRt(double rt) { rt_ = rt; }

  double getBaseInte() const { return base_inte_;}
  void setBaseInte(double base_inte) {base_inte_ = base_inte;}

 private:
  int row_id_;
  int raw_spec_id_;
  int scan_num_;
  double rt_;
  double base_inte_;
};

typedef std::shared_ptr<MsMapRowHeader> MsMapRowHeaderPtr;
typedef std::vector<MsMapRowHeaderPtr> MsMapRowHeaderPtrVec;

}

#endif 
