
#ifndef PROT_NEUTRAL_LOSS_HPP_
#define PROT_NEUTRAL_LOSS_HPP_

#include <string>
#include <vector>
#include <memory>

namespace prot {

#define NEUTRAL_LOSS_NONE "NONE"

class NeutralLoss {

 public:
  NeutralLoss(const std::string &name, double mass);
  std::string getName(){return name_;}
  double getMass(){return mass_;}

 private:
  std::string name_;
  double mass_;
};

typedef std::shared_ptr<NeutralLoss> NeutralLossPtr;
typedef std::vector<NeutralLossPtr> NeutralLossPtrVec;

/* neutral loss factory */
class NeutralLossFactory {
 public:
  static void initFactory(const std::string &file_name);
  static NeutralLossPtrVec getBaseNeutralLossPtrVec() {
    return neutral_loss_ptr_vec_;}
  static NeutralLossPtr getBaseNeutralLossPtrByName(const std::string &name);

  static NeutralLossPtr getNeutralLossPtr_NONE () {
    return getBaseNeutralLossPtrByName(NEUTRAL_LOSS_NONE);
  }

 private:
  static NeutralLossPtrVec neutral_loss_ptr_vec_;
};

} /* namespace prot */

#endif /* NEUTRALLOSS_HPP_ */
