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

#include "prsm/prsm_sql_writer.hpp"

#include <sqlite3.h>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "common/util/str_util.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "seq/alter_type.hpp"
#include "seq/fasta_seq.hpp"
#include "seq/mass_shift.hpp"
#include "seq/proteoform.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace {

void bindText(sqlite3_stmt* stmt, int idx, const std::string& text) {
  sqlite3_bind_text(stmt, idx, text.c_str(), -1, SQLITE_TRANSIENT);
}

// FDR values are -1 when they were not computed; store those as NULL.
void bindFdr(sqlite3_stmt* stmt, int idx, double fdr) {
  if (fdr >= 0) {
    sqlite3_bind_double(stmt, idx, fdr);
  } else {
    sqlite3_bind_null(stmt, idx);
  }
}

}  // namespace

PrsmSqlWriter::PrsmSqlWriter(const PrsmParaPtr& prsm_para_ptr, sqlite3* sql_db)
    : prsm_para_ptr_(prsm_para_ptr), sql_db_(sql_db) {
  std::string db_file_name = prsm_para_ptr_->getSearchDbFileNameWithFolder();
  search_match_ptr_ = std::make_shared<SearchFastaMatch>(db_file_name);
}

void PrsmSqlWriter::createPrsmTable(const std::string& table_name) {
  sql_util::execSql(sql_db_, "DROP TABLE IF EXISTS " + table_name + ";");
  sql_util::execSql(sql_db_, "CREATE TABLE " + table_name +
                                 "("
                                 "prsm_id INTEGER PRIMARY KEY,"
                                 "spectrum_id INT NOT NULL,"
                                 "fragmentation TEXT NOT NULL,"
                                 "scans TEXT NOT NULL,"
                                 "retention_time REAL NOT NULL,"  // seconds
                                 "mass_num INT NOT NULL,"
                                 "charge INT NOT NULL,"
                                 "precursor_mass REAL NOT NULL,"
                                 "adjusted_precursor_mass REAL NOT NULL,"
                                 "proteoform_id INT NOT NULL,"
                                 "proteoform_intensity REAL NOT NULL,"
                                 "feature_id INT,"
                                 "feature_intensity REAL,"
                                 "feature_score REAL,"
                                 "feature_apex_time REAL,"  // seconds
                                 "protein_hit_num INT NOT NULL,"
                                 "protein_name TEXT NOT NULL,"
                                 "protein_description TEXT NOT NULL,"
                                 "first_residue INT NOT NULL,"
                                 "last_residue INT NOT NULL,"
                                 "special_amino_acids TEXT NOT NULL,"
                                 "protein_sequence TEXT NOT NULL,"
                                 "prev_amino_acid TEXT NOT NULL,"
                                 "proteoform TEXT NOT NULL,"
                                 "next_amino_acid TEXT NOT NULL,"
                                 "proteoform_mass REAL NOT NULL,"
                                 "n_terminal_form TEXT NOT NULL,"
                                 "fixed_ptms TEXT NOT NULL,"
                                 "unexpected_mod_num INT NOT NULL,"
                                 "unexpected_mods TEXT NOT NULL,"
                                 "variable_ptm_num INT NOT NULL,"
                                 "variable_ptms TEXT NOT NULL,"
                                 "mi_score TEXT NOT NULL,"
                                 "matched_peak_num INT NOT NULL,"
                                 "matched_fragment_num INT NOT NULL,"
                                 "e_value REAL NOT NULL,"
                                 "spectrum_fdr REAL,"
                                 "proteoform_fdr REAL,"
                                 "protein_fdr REAL);");
}

void PrsmSqlWriter::write(const std::string& input_file_ext,
                          const std::string& table_name, bool write_details) {
  PrsmPtrVec prsm_list =
      prsm_reader_util::readPrsmsWithSpectra(prsm_para_ptr_, input_file_ext);

  sql_util::execSql(sql_db_, "BEGIN;");
  createPrsmTable(table_name);
  sqlite3_stmt* prsm_stmt = sql_util::prepareSql(
      sql_db_, "INSERT INTO " + table_name +
                   " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
                   "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
                   "?, ?, ?, ?);");
  sqlite3_stmt* shift_stmt = nullptr;
  sqlite3_stmt* match_stmt = nullptr;
  if (write_details) {
    sql_util::execSql(sql_db_, "DROP TABLE IF EXISTS prsm_mass_shift;");
    sql_util::execSql(sql_db_,
                      "CREATE TABLE prsm_mass_shift("
                      "prsm_id INT NOT NULL,"
                      "shift_index INT NOT NULL,"
                      "left_position INT NOT NULL,"   // 0-based residue index
                      "right_position INT NOT NULL,"  // in the proteoform
                      "mass REAL NOT NULL,"
                      "type TEXT NOT NULL,"
                      "annotation TEXT NOT NULL,"
                      "PRIMARY KEY (prsm_id, shift_index));");
    shift_stmt = sql_util::prepareSql(
        sql_db_, "INSERT INTO prsm_mass_shift VALUES (?, ?, ?, ?, ?, ?, ?);");
    sql_util::execSql(sql_db_, "DROP TABLE IF EXISTS prsm_protein_match;");
    sql_util::execSql(sql_db_,
                      "CREATE TABLE prsm_protein_match("
                      "prsm_id INT NOT NULL,"
                      "protein_name TEXT NOT NULL,"
                      "protein_description TEXT NOT NULL,"
                      "first_residue INT NOT NULL,"
                      "last_residue INT NOT NULL,"
                      "PRIMARY KEY (prsm_id, protein_name));");
    match_stmt = sql_util::prepareSql(
        sql_db_, "INSERT INTO prsm_protein_match VALUES (?, ?, ?, ?, ?);");
  }

  for (const PrsmPtr& prsm_ptr : prsm_list) {
    writePrsm(prsm_stmt, prsm_ptr);
    if (write_details) {
      writeMassShifts(shift_stmt, prsm_ptr);
      writeProteinMatches(match_stmt, prsm_ptr);
    }
  }
  sql_util::execSql(sql_db_, "COMMIT;");
  sqlite3_finalize(prsm_stmt);
  sqlite3_finalize(shift_stmt);
  sqlite3_finalize(match_stmt);
}

void PrsmSqlWriter::writePrsm(sqlite3_stmt* stmt, const PrsmPtr& prsm_ptr) {
  // The PrSM may cover a group of spectra; list their ids/activations/scans.
  std::string activations;
  std::string scans;
  int mass_num = 0;
  const DeconvMsPtrVec& ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
    MsHeaderPtr header_ptr = ms_ptr_vec[i]->getMsHeaderPtr();
    activations += header_ptr->getActivationPtr()->getName() + " ";
    scans += header_ptr->getScansString() + " ";
    mass_num += ms_ptr_vec[i]->size();
  }
  str_util::trim(activations);
  str_util::trim(scans);
  MsHeaderPtr first_header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();

  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  int start_pos = form_ptr->getStartPos();
  int end_pos = form_ptr->getEndPos();
  std::vector<std::pair<FastaSeqPtr, int>> matches =
      search_match_ptr_->process(prsm_ptr);

  int idx = 1;
  sqlite3_bind_int(stmt, idx++, prsm_ptr->getPrsmId());
  sqlite3_bind_int(stmt, idx++, prsm_ptr->getSpectrumId());
  bindText(stmt, idx++, activations);
  bindText(stmt, idx++, scans);
  sqlite3_bind_double(stmt, idx++, first_header_ptr->getRetentionTime());
  sqlite3_bind_int(stmt, idx++, mass_num);
  sqlite3_bind_int(stmt, idx++, first_header_ptr->getFirstPrecCharge());
  sqlite3_bind_double(stmt, idx++, prsm_ptr->getOriPrecMass());
  sqlite3_bind_double(stmt, idx++, prsm_ptr->getAdjustedPrecMass());
  sqlite3_bind_int(stmt, idx++, form_ptr->getProteoClusterId());
  sqlite3_bind_double(stmt, idx++, form_ptr->getProteoInte());
  if (prsm_ptr->getFracFeatureInte() > 0) {
    sqlite3_bind_int(stmt, idx++, prsm_ptr->getFracFeatureId());
    sqlite3_bind_double(stmt, idx++, prsm_ptr->getFracFeatureInte());
    sqlite3_bind_double(stmt, idx++, prsm_ptr->getFracFeatureScore());
    sqlite3_bind_double(stmt, idx++, prsm_ptr->getFracFeatureApexTime());
  } else {
    for (int i = 0; i < 4; i++) {
      sqlite3_bind_null(stmt, idx++);
    }
  }
  sqlite3_bind_int(stmt, idx++, static_cast<int>(matches.size()));
  bindText(stmt, idx++, form_ptr->getSeqName());
  bindText(stmt, idx++, form_ptr->getSeqDesc());
  sqlite3_bind_int(stmt, idx++, start_pos + 1);
  sqlite3_bind_int(stmt, idx++, end_pos + 1);
  bindText(stmt, idx++,
           form_ptr->getFastaSeqPtr()->getAcidReplaceStr(start_pos, end_pos));
  bindText(stmt, idx++,
           form_ptr->getFastaSeqPtr()->getSubSeq(start_pos, end_pos));
  bindText(stmt, idx++, form_ptr->getPrevAminoAcid());
  bindText(stmt, idx++, form_ptr->getProteoformMatchSeq());
  bindText(stmt, idx++, form_ptr->getNextAminoAcid());
  sqlite3_bind_double(stmt, idx++, form_ptr->getMass());
  bindText(stmt, idx++, form_ptr->getProtModPtr()->getType());
  bindText(stmt, idx++, form_ptr->getAlterStr(AlterType::FIXED));
  sqlite3_bind_int(stmt, idx++, form_ptr->getAlterNum(AlterType::UNEXPECTED));
  bindText(stmt, idx++, form_ptr->getAlterStr(AlterType::UNEXPECTED));
  sqlite3_bind_int(stmt, idx++, form_ptr->getAlterNum(AlterType::VARIABLE));
  bindText(stmt, idx++, form_ptr->getAlterStr(AlterType::VARIABLE));
  bindText(stmt, idx++, form_ptr->getMIScore());
  sqlite3_bind_int(stmt, idx++, static_cast<int>(prsm_ptr->getMatchPeakNum()));
  sqlite3_bind_int(stmt, idx++, static_cast<int>(prsm_ptr->getMatchFragNum()));
  sqlite3_bind_double(stmt, idx++, prsm_ptr->getEValue());
  bindFdr(stmt, idx++, prsm_ptr->getFdr());
  bindFdr(stmt, idx++, prsm_ptr->getProteoformFdr());
  bindFdr(stmt, idx++, prsm_ptr->getProteinFdr());
  sql_util::stepAndReset(sql_db_, stmt);
}

void PrsmSqlWriter::writeMassShifts(sqlite3_stmt* stmt,
                                    const PrsmPtr& prsm_ptr) {
  const MassShiftPtrVec& shifts =
      prsm_ptr->getProteoformPtr()->getMassShiftPtrVec();
  for (size_t i = 0; i < shifts.size(); i++) {
    sqlite3_bind_int(stmt, 1, prsm_ptr->getPrsmId());
    sqlite3_bind_int(stmt, 2, static_cast<int>(i));
    sqlite3_bind_int(stmt, 3, shifts[i]->getLeftBpPos());
    sqlite3_bind_int(stmt, 4, shifts[i]->getRightBpPos());
    sqlite3_bind_double(stmt, 5, shifts[i]->getMassShift());
    bindText(stmt, 6, shifts[i]->getTypePtr()->getName());
    bindText(stmt, 7, shifts[i]->getAnnoStr());
    sql_util::stepAndReset(sql_db_, stmt);
  }
}

void PrsmSqlWriter::writeProteinMatches(sqlite3_stmt* stmt,
                                        const PrsmPtr& prsm_ptr) {
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  std::vector<std::pair<FastaSeqPtr, int>> matches =
      search_match_ptr_->process(prsm_ptr);
  for (const auto& match : matches) {
    const FastaSeqPtr& seq_ptr = match.first;
    if (seq_ptr->getName() == form_ptr->getSeqName()) {
      continue;
    }
    int seq_pos = match.second;
    sqlite3_bind_int(stmt, 1, prsm_ptr->getPrsmId());
    bindText(stmt, 2, seq_ptr->getName());
    bindText(stmt, 3, seq_ptr->getDesc());
    sqlite3_bind_int(stmt, 4, seq_pos + 1);
    sqlite3_bind_int(stmt, 5, seq_pos + form_ptr->getLen());
    sql_util::stepAndReset(sql_db_, stmt);
  }
}

}  // namespace toppic
