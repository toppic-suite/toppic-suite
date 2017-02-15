// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_SPEC_MS_HPP_
#define PROT_SPEC_MS_HPP_

#include <sstream>

#include "spec/ms_header.hpp"

namespace prot {

template <class T>
class Ms {
 public:
  Ms() {};

  Ms(MsHeaderPtr header_ptr) {header_ptr_ = header_ptr;}

  Ms(MsHeaderPtr header_ptr, const std::vector<T> &peak_ptr_list) {
    header_ptr_ = header_ptr;
    peak_ptr_list_ = peak_ptr_list;
  }

  /**
   * Removes precursor mass. In ETD data, MSMS may contain a high precursor
   * mass peak. So we use the following to remove it.
   */
  void rmPrec(double tolerance) {
    peak_ptr_list_ = rmPeaks(peak_ptr_list_, header_ptr_->getPrecSpMz(), 
                             tolerance);
  }

  void recalibrate(double recal) {
    for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
      double new_mass = (1 + recal) * peak_ptr_list_[i]->getPosition();
      peak_ptr_list_[i]->setPosition(new_mass);
    }
  }

  std::string toString() {
    std::string header_str = header_ptr_->toString();
    std::stringstream tmp;
    for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
      tmp << i << " " << peak_ptr_list_[i]->getPosition() 
          << " " << peak_ptr_list_[i]->getIntensity() << "\n";
    }
    return header_str + tmp.str();
  }

  MsHeaderPtr getMsHeaderPtr() {return header_ptr_;}

  void setHeaderPtr(MsHeaderPtr header_ptr) {header_ptr_ = header_ptr;}

  size_t size() {return peak_ptr_list_.size();}

  T getPeakPtr(int i) {return peak_ptr_list_[i];}
  
  std::vector<T> getPeakPtrVec() {return peak_ptr_list_;}

 private:
  MsHeaderPtr header_ptr_;
  std::vector<T> peak_ptr_list_;
};

}
#endif
