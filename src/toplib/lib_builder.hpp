// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef TOPPIC_TOPLIB_LIB_BUILDER_HPP_
#define TOPPIC_TOPLIB_LIB_BUILDER_HPP_

#include <string>
#include <vector>

#include "toplib/lib_spectrum.hpp"

namespace toppic {

namespace toplib {

// How a cluster's representative spectrum is made: the merged spectrum of
// all members ("average") or the member with the lowest E-value ("single").
enum class LibRepMethod { kAverage, kSingle };

// The number of peaks a representative keeps.
inline constexpr int kLibPeakNum = 50;

// One representative spectrum of the library, with its decoy.
struct LibRep {
  double precursor_mass = 0;
  double precursor_intensity = 0;
  int precursor_charge = 0;
  // at most kLibPeakNum peaks by increasing mass
  std::vector<double> masses;
  std::vector<double> intensities;
  std::vector<int> charges;
  // log2 of the intensities, scaled to unit length
  std::vector<double> norm_intensities;
  // the decoy's masses (in the order of masses)
  std::vector<double> decoy_masses;
  // one value per LibLibrary::ptm_columns: the identification and
  // spectrum information of the cluster's best member
  std::vector<LibValue> ptm_values;
  // false when the cluster has no identification or loses it to a better
  // cluster of the same proteoform and precursor charge; such a
  // representative is not written
  bool identified = false;
};

// The library before it is written: one representative per cluster, in
// cluster order (representative_id = position + 1).
struct LibLibrary {
  std::vector<LibIdColumn> ptm_columns;
  std::vector<LibRep> reps;
  // the members of each cluster (positions in LibData::spectra, increasing)
  std::vector<std::vector<int>> cluster_members;
  // counts for the log
  int input_spectrum_num = 0;
  int valid_spectrum_num = 0;
  int precursor_cluster_num = 0;
  int unidentified_cluster_num = 0;
  int inconsistent_cluster_num = 0;
};

namespace lib_builder {

// ms_first_intensity_remove: keeps the spectra whose first precursor has a
// nonzero intensity that is at least cutoff of the summed intensity of all
// the spectrum's precursors (the precursor-to-precursor intensity ratio).
// A cutoff of 0 keeps every spectrum, as the Python script skips the
// filter then.
void filterPpir(LibData& data, double cutoff);

// ms_rep_library_building: removes the spectra with fewer than two peaks or
// a zero precursor intensity, clusters the rest by precursor charge,
// precursor mass (complete linkage, 2.2 Da) and the cosine similarity of
// their 50 most intense peaks (average linkage, distance 0.7), makes one
// representative per cluster, gives it the identification of its member
// with the lowest E-value (keeping only the best cluster per proteoform and
// precursor charge) and a decoy spectrum. Spectra identified as different
// proteoforms (proteoform id and protein accession) are never clustered
// together, so every proteoform keeps a representative; python_compatible
// lets them be and keeps one cluster per proteoform id and charge, as the
// Python script does.
LibLibrary build(const LibData& data, LibRepMethod method,
                 bool python_compatible = false);

}  // namespace lib_builder

namespace lib_sql_writer {

// Writes the representatives with an identification, targets then decoys,
// into the table target_decoy_spectra_representatives of a SQLite database
// (the table is replaced, other tables are kept), with the columns and
// column types the Python script's pandas data frame gets.
void write(const LibLibrary& lib, const std::string& db_file_name);

}  // namespace lib_sql_writer

}  // namespace toplib

}  // namespace toppic

#endif
