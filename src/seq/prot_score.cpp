//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "seq/prot_score.hpp"

namespace toppic {

ProtScore::ProtScore(int id, int score, 
                     double n_term_shift,
                     double c_term_shift):  
  id_(id),
  score_(score),
  n_term_shift_(n_term_shift),
  c_term_shift_(c_term_shift) {}

} /* namespace toppic */
