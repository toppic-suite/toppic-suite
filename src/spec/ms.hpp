#ifndef PROT_MS_HPP_
#define PROT_MS_HPP_

#include "ms_header.hpp"

namespace prot {

template <class T>
class Ms {
 public:
  Ms() {};

  Ms(MsHeaderPtr header_ptr) {
    header_ptr_ = header_ptr;
  }

  Ms(MsHeaderPtr header_ptr, std::vector<T> peak_ptr_list) {
    header_ptr_ = header_ptr;
    peak_ptr_list_ = peak_ptr_list;
  }

	/**
	 * Removes precursor mass. In ETD data, MSMS may contain a high precursor
	 * mass peak. So we use the following to remove it.
	 */
  void rmPrec(double tolerance);

	void recalibrate(double recal);

  std::string toString();

	MsHeaderPtr getHeaderPtr() {return header_ptr_;}

	void setHeaderPtr(MsHeaderPtr header_ptr) {header_ptr = header_ptr_;}

  unsigned int size() {return peak_ptr_list_.size();}

  T getPeakPtr(int i) {return peak_ptr_list_[i];}

 private:
  MsHeaderPtr header_ptr_;
  std::vector<T> peak_ptr_list_;
};

}
#endif
