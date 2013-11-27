#include<string>
#include<sstream>

#include "ms.hpp"

namespace prot {

template <class T>
Ms<T>::Ms(MsHeaderPtr header_ptr) {
  header_ptr_ = header_ptr;
}

template <class T>
Ms<T>::Ms(MsHeaderPtr header_ptr, std::vector<T> peak_ptr_list) {
  header_ptr_ = header_ptr;
  peak_ptr_list_ = peak_ptr_list;
}

template <class T>
void Ms<T>::rmPrec(double tolerance) {
  peak_ptr_list_ = rmPeaks(peak_ptr_list_, header_ptr_->getPrecSpMz(), tolerance);
}

template <class T>
void Ms<T>::recalibrate(double recal) {
  for (unsigned int i = 0; i < peak_ptr_list_.size(); i++) {
    double new_mass = (1 + recal) * peak_ptr_list_[i]->getPosition();
    peak_ptr_list_[i]->setPosition(new_mass);
  }
}

template <class T>
std::string Ms<T>::toString() {
  std::string header_str = header_ptr_->toString();
  std::stringstream tmp;
  for (unsigned int i = 0; i < peak_ptr_list_.size(); i++) {
    tmp << i << " " << peak_ptr_list_[i]->getPosition() << " " << peak_ptr_list_[i]->getIntensity() << "\n";
  }
  return header_str + tmp.str();
}

}
