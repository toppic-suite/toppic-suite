#ifndef PROT_SUMMARY_HPP_
#define PROT_SUMMARY_HPP_

#include "boost/date_time/posix_time/posix_time.hpp"

#include "console/argument.hpp"

namespace prot {

class Summary {
 public:
  static void outputSummary(Argument arguments, 
                            boost::posix_time::ptime start_time,
                            boost::posix_time::ptime stop_time);
};
}

#endif /* ARGUMENT_HPP_ */
