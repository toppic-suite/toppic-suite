
#ifndef PROT_ION_TYPE_H_
#define PROT_ION_TYPE_H_

#include <string>
#include <vector>
#include <memory>

namespace prot {

class IonType;
typedef std::shared_ptr<IonType> IonTypePtr;

class IonType {
 public: 
  static const IonTypePtr B;
  static const IonTypePtr Y;
  static const IonTypePtr C;
  static const IonTypePtr Z_DOT;

  IonType(std::string name, bool n_term, double shift);

  std::string getName() {return name_;}

  bool isNTerm() {return n_term_;}

  double getShift() {return shift_;}

	/** Gets IonType by name */
	static IonTypePtr getIonTypePtrByName(std::string name);

 private:
  /** ion name */
  std::string name_;
  /** A B C are n-terminal ions and X Y Z are non-n-terminal ions */
  bool n_term_;
  /**
   * Shift stands for the shift of the ion compared to M. For example, the
   * shift for B ion is 0, and the shift for Y ion is 18 (chrg 0);
   */
  double shift_;
};


}

#endif
