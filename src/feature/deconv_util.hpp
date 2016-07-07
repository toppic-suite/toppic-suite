#ifndef PROT_FEATURE_DECONV_UTIL_HPP_
#define PROT_FEATURE_DECONV_UTIL_HPP_

#include <memory>
#include <vector>

namespace prot {

class IntvDens {
 public:
  IntvDens(double bgn, double end, int num, double perc): 
      bgn_(bgn), 
      end_(end),
      num_(num),
      perc_(perc) {
      }

  double getBgn() {return bgn_;}   
  double getEnd() {return end_;}   
  double getMiddle() {return (bgn_ + end_) / 2;}   
  int getNum() {return num_;}   
  double getPerc() {return perc_;}   

  private:
   double bgn_, end_;
   int num_;
   double perc_;
};

typedef std::shared_ptr<IntvDens> IntvDensPtr;
typedef std::vector<IntvDensPtr> IntvDensPtrVec;

class DeconvUtil {
 public:
  static IntvDensPtrVec getDensity(std::vector<double> &inte);

  static void outputDens(IntvDensPtrVec &dens);

  static int getMaxPos(IntvDensPtrVec &dens);

  static double getBaseLine(std::vector<double> &inte);

 private:
  static double intv_width_;
};


}
#endif
