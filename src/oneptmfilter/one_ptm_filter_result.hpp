#ifndef ONE_PTM_FILTER_RESULT_HPP_
#define ONE_PTM_FILTER_RESULT_HPP_

#include "base/proteoform.hpp"
#include "base/base_data.hpp"

namespace prot {

class OnePtmFilterResult {

 private:
  int proteo_id_;
  int score_;
  int pref_pos_;
  int suff_pos_;
}

}

#endif /* ONE_PTM_FILTER_RESULT_HPP_ */
