#ifndef PROT_FEATURE_DECONV_PARA_HPP_
#define PROT_FEATURE_DECONV_PARA_HPP_

#include <string>
#include <memory>

namespace prot {

enum OutputType {OUTPUT_MGF, OUTPUT_MSALIGN};
enum InputType {INPUT_MGF, INPUT_MZXML};

class DeconvPara {
 public:

  void setDataFileName(std::string &file_name) {data_file_name_ == file_name;}

  std::string getDataFileName() {return data_file_name_;}

  int setOutputType (std::string &format);

  int setInputType (std::string &format);

  std::string data_file_name_;
  InputType input_type_;
  OutputType output_type_;

  bool refine_prec_mass_;
  int ms_level_; //ms level
  bool missing_level_one_;

  int max_charge_;
  double max_mass_;
  double tolerance;
  double sn_ratio_;
  bool keep_unused_peaks_;

  bool output_multiple_mass_ = false; 

  double prec_window_;
};

typedef std::shared_ptr<DeconvPara> DeconvParaPtr;

}
#endif
