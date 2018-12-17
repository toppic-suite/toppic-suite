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


#ifndef PROT_FEATURE_DECONV_UTIL_HPP_
#define PROT_FEATURE_DECONV_UTIL_HPP_

#include <memory>
#include <vector>

namespace toppic {

namespace deconv_util {

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

IntvDensPtrVec getDensity(const std::vector<double> &inte);

void outputDens(const IntvDensPtrVec &dens);

int getMaxPos(const IntvDensPtrVec &dens);

double getBaseLine(const std::vector<double> &inte);

}  // namespace deconv_util 

}  // namespace toppic
#endif
