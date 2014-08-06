#ifndef PROT_MS_HPP_
#define PROT_MS_HPP_

#include <sstream>

#include "spec/ms_header.hpp"

namespace prot {

template <class T>
class Ms {
 public:
  Ms() {};

  Ms(MsHeaderPtr header_ptr) {
    header_ptr_ = header_ptr;
  }

  Ms(MsHeaderPtr header_ptr, const std::vector<T> &peak_ptr_list) {
    header_ptr_ = header_ptr;
    peak_ptr_list_ = peak_ptr_list;
  }

  Ms(MsHeaderPtr header_ptr, const std::vector<T> &peak_ptr_list, double ppo) {
    header_ptr_ = header_ptr;
    peak_ptr_list_ = peak_ptr_list;
    header_ptr_->setErrorToleranceByPpo(ppo);
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

  MsHeaderPtr getHeaderPtr() {return header_ptr_;}

  void setHeaderPtr(MsHeaderPtr header_ptr) {header_ptr = header_ptr_;}

  size_t size() {return peak_ptr_list_.size();}

  T getPeakPtr(int i) {return peak_ptr_list_[i];}
  
  std::vector<T> getPeakPtrVec() {return peak_ptr_list_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
        xercesc::DOMElement* element = xml_doc->createElement("ms");
        header_ptr_->appendXml(xml_doc,element);
        parent->appendChild(element);
  }

 private:
  MsHeaderPtr header_ptr_;
  std::vector<T> peak_ptr_list_;
};

}
#endif
