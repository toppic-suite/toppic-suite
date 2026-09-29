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

#ifndef TOPPIC_TOPLIB_LIB_MSALIGN_READER_HPP_
#define TOPPIC_TOPLIB_LIB_MSALIGN_READER_HPP_

#include <string>

#include "toplib/lib_spectrum.hpp"

namespace toppic {

namespace toplib {

namespace lib_msalign_reader {

// A TopFD msalign file and a filtered TopPIC PrSM table, as TopLib's
// ms_library_building.py reads them (msalign_file.read_msalign and
// get_spectra_by_file): header fields are found by substring, spectra
// without masses or with a zero precursor mass are skipped, masses above the
// precursor mass are dropped, and the table rows are joined by spectrum id.
// The table's header is its first line with a Spectrum ID column (TopLib
// wants it on the first line; TopPIC writes its parameters above it).
LibData readMsalignAndPrsms(const std::string& msalign_file_name,
                            const std::string& prsm_file_name);

// A library msalign file (v1.2.1: the identification in the spectrum
// headers, e.g. PROTEIN_ACCESSION, PROTEOFORM, E_VALUE; MS2_ID, MS2_SCAN,
// PRECURSOR_MONOISOTOPIC_MASS, ...). The fields are mapped to the columns
// the msalign file and the PrSM table give; the proteoform id is the index
// of the spectrum's distinct (accession, proteoform) pair in file order and
// the last residue is the first residue plus the proteoform's residue count
// minus one.
LibData readLibraryMsalign(const std::string& msalign_file_name);

}  // namespace lib_msalign_reader

}  // namespace toplib

}  // namespace toppic

#endif
