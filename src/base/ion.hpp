
#ifndef PROT_ION_HPP_
#define PROT_ION_HPP_

#include <string>

#include "base/ion_type.hpp"
#include "base/neutral_loss.hpp"

namespace prot {

class Ion {
 public:
  Ion(int charge, int pos, int display_pos, 
      IonTypePtr ion_type_ptr, 
      NeutralLossPtr neutral_loss_ptr);

  int getCharge() {return charge_;}
  int getPos() {return pos_;}
  int getDisplayPos() {return display_pos_;}

  IonTypePtr getIonTypePtr() {return ion_type_ptr_;}

  std::string getDisplayName(){
    std::string pos_str = std::to_string(display_pos_);
    return ion_type_ptr_ -> getName() + pos_str;
  }

 private:
  int charge_;
  int pos_;
  int display_pos_;
  IonTypePtr ion_type_ptr_;
  NeutralLossPtr neutral_loss_ptr_;
};

typedef std::shared_ptr<Ion> IonPtr;
typedef std::vector<IonPtr> IonPtrVec;

} /* namespace prot */

#endif /* ION_HPP_ */
