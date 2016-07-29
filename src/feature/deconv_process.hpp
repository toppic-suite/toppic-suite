#ifndef PROT_FEATURE_DECONV_PROCESS_HPP_
#define PROT_FEATURE_DECONV_PROCESS_HPP_

#include <string>
#include "feature/deconv_para.hpp"
#include "feature/deconv_one_sp.hpp"
#include "feature/feature_mng.hpp"
#include "feature/feature_ms_reader.hpp"

namespace prot {

class DeconvProcess {
 public:
  DeconvProcess(DeconvParaPtr para_ptr) {para_ptr_ = para_ptr;}
	int getResult() {return result_;}
  std::string getMsg() {return msg_;}

  void process();

  void processSp(DeconvOneSpPtr deconv_ptr, FeatureMsReaderPtr reader_ptr, 
                 std::ofstream &os1, std::ofstream &os2);

 private:
  DeconvParaPtr para_ptr_;
  std::string msg_ = "";
  int result_ = 0;

  void copyParameters(FeatureMngPtr mng_ptr);
  void printParameter(FeatureMngPtr mng_ptr);

  void updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num);
};

}
#endif
