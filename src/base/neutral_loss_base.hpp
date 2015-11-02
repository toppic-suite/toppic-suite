#ifndef PROT_BASE_NEUTRAL_LOSS_BASE_HPP_
#define PROT_BASE_NEUTRAL_LOSS_BASE_HPP_

#include <string>
#include <vector>
#include <memory>

#include "base/neutral_loss.hpp"

namespace prot {

class NeutralLossBase {
 public:
  static void initBase(const std::string &file_name);
  static NeutralLossPtrVec getBaseNeutralLossPtrVec() {
    return neutral_loss_ptr_vec_;}

  static NeutralLossPtr getNeutralLossPtrByName(const std::string &name);

  static NeutralLossPtr getNeutralLossPtr_NONE () {
    return neutral_loss_ptr_NONE_;
  }

  static std::string getName_NONE() {return "NONE";}

 private:
  static NeutralLossPtrVec neutral_loss_ptr_vec_;
  static NeutralLossPtr neutral_loss_ptr_NONE_;
};

} /* namespace prot */

#endif /* NEUTRALLOSS_HPP_ */
