/*
 * ion.cpp
 *
 *  Created on: Nov 25, 2013
 *      Author: xunlikun
 */

#include "base/ion.hpp"

namespace prot {

Ion::Ion(int charge,int pos,int display_pos,
         IonTypePtr ion_type_ptr,
         NeutralLossPtr neutral_loss_ptr){
  charge_ = charge;
  pos_=pos;
  display_pos_=display_pos;
  ion_type_ptr_=ion_type_ptr;
  neutral_loss_ptr_=neutral_loss_ptr;
}

} /* namespace prot */
