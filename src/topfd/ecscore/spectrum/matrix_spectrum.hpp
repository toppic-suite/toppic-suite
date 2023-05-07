//
// Created by abbash on 8/22/22.
//

#ifndef TOPPIC_TOPFD_ECSCORE_MATRIX_SPECTRUM_HPP
#define TOPPIC_TOPFD_ECSCORE_MATRIX_SPECTRUM_HPP

#include <memory>
#include <vector>

namespace toppic {

class MatrixSpectrum {
 public:
  MatrixSpectrum(int spec_id, int scan_num, double rt); 

  int getSpecId() const { return spec_id_; }
  void setSpecId(int specId) { spec_id_ = specId; }

  int getScanNum() const { return scan_num_; }
  void setScanNum(int scanNum) { scan_num_ = scanNum; }

  double getRt() const { return rt_; }
  void setRt(double rt) { rt_ = rt; }

 private:
  int spec_id_;
  int scan_num_;
  double rt_;

};

typedef std::shared_ptr<MatrixSpectrum> MatrixSpectrumPtr;
typedef std::vector<MatrixSpectrumPtr> MatrixSpectrumPtrVec;

}

#endif //TOPPIC_SPECTRUM_HPP
