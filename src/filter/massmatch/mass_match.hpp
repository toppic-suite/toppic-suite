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

#ifndef TOPPIC_FILTER_MASS_MATCH_MATCH_MATCH_HPP_
#define TOPPIC_FILTER_MASS_MATCH_MATCH_MATCH_HPP_

#include <memory>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

namespace toppic {

class MassMatch {
 public:
  MassMatch(){};
  MassMatch(std::vector<std::vector<int>> &mass_2d,
            std::vector<std::vector<double>> &real_shift_2d,
            std::vector<std::vector<int>> &pos_2d,
            double max_proteoform_mass, double scale);
            MassMatch(std::vector<std::vector<int>> &mass_2d,
            std::vector<std::vector<double>> &real_shift_2d,
            std::vector<std::vector<int>> &pos_2d,
            double max_proteoform_mass, double scale, bool prm);

  void compScores(const std::vector<std::pair<int, int>> &pref_mass_errors,
                  std::vector<short> &scores);

  void compScores(const std::vector<std::pair<int, int>> &pref_mass_errors,
                  int start, double shift, std::vector<short> &scores);

  void compMatchScores(const std::vector<std::pair<int, int>> &pref_mass_errors,
                       const std::pair<int, int> &prec_minus_water_mass_error,
                       std::vector<short> &scores);

  //void serializeMassMatch();
  void serializeMassMatch(std::string fName, std::string block_str, std::string dirName);
  //void deserializeMassMatch(MassMatch **m);
  void deserializeMassMatch(MassMatch **m);

  int getProtNum() {return proteo_num_;}

  int getRowNum() {return row_num_;}

  int getColNum() {return col_num_;}

  bool getPrm() {return prm_;}
  
  std::string getFileName(){return file_name_;}

  std::string getDirName(){return dir_name_;}
  
  static int getPrecursorMatchScore() {return 10000;}

  const std::vector<int>& getProteoRowBegins() {return proteo_row_begins_;}

  const std::vector<int>& getProteoRowEnds() {return proteo_row_ends_;}

  const std::vector<double>& getTruncShifts() {return trunc_shifts_;}
  
  //set file name
  void setfileName(std::string name){file_name_ = name;}
  void setDirName(std::string dir){dir_name_ = dir;}

 private:
 //for serialization
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
      ar & scale_;
      ar & proteo_num_;
      ar & col_num_; 
      ar & row_num_;
    
      ar & proteo_row_begins_;
      ar & proteo_row_ends_;
      ar & row_proteo_ids_;
      ar & trunc_shifts_;

      ar & col_index_begins_;
      ar & col_index_ends_;
      ar & col_indexes_;
      
  }

  double scale_;

  int proteo_num_;

  int col_num_;
  int row_num_;
  bool prm_;

  std::string file_name_;
  std::string dir_name_;

  // the first row of each proteoform
  std::vector<int> proteo_row_begins_;
  std::vector<int> proteo_row_ends_;

  // the proteoform id of each row
  std::vector<int> row_proteo_ids_;
  std::vector<double> trunc_shifts_;

  std::vector<int> col_index_begins_;
  std::vector<int> col_index_ends_;
  std::vector<int> col_indexes_;

  void initProteoformBeginEnds(std::vector<std::vector<double>> &shift_2d);

  void initIndexes(std::vector<std::vector<int>> &mass_2d,
                   std::vector<std::vector<double>> &real_shift_2d,
                   std::vector<std::vector<int>> &pos_2d);

  void compColumnMatchNums(std::vector<std::vector<int>> &mass_2d,
                           std::vector<std::vector<int>> &shift_2d,     
                           std::vector<std::vector<int>> &pos_2d,
                           std::vector<int> &col_match_nums);

  void fillColumnIndex(std::vector<std::vector<int>> &mass_2d,
                       std::vector<std::vector<int>> &shift_2d,
                       std::vector<std::vector<int>> &pos_2d,
                       std::vector<int> &col_index_pnts);
};
    typedef std::shared_ptr<MassMatch> MassMatchPtr;

} /* namespace toppic */

#endif
