/*
 * diagonal_header.hpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#ifndef DIAGONAL_HEADER_HPP_
#define DIAGONAL_HEADER_HPP_

#include "base/trunc.hpp"
#include "base/prot_mod.hpp"
#include "base/ptm.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"

namespace prot {
class DiagonalHeader;
typedef std::shared_ptr<DiagonalHeader> DiagonalHeaderPtr;
typedef std::vector<DiagonalHeaderPtr> DiagonalHeaderPtrVec;
typedef std::vector<DiagonalHeaderPtrVec> DiagonalHeaderPtrVec2D;
typedef std::vector<DiagonalHeaderPtrVec2D> DiagonalHeaderPtrVec3D;

class DiagonalHeader {
public:
	DiagonalHeader(double n_term_shift,bool n_strict,bool c_strict,bool n_trunc,bool c_trunc);
	DiagonalHeaderPtr clone();
	void changeNTermShift(double s){
		prot_N_term_shift_ -= s;
		pep_N_term_shift_ -= s;
		prot_C_term_shift_ +=s;
		pep_C_term_shift_+=s;
	}

	int getTruncFirstTesPos(){return trunc_first_res_pos_;}

	int getMatchFirstResPos() const {
		return match_first_res_pos_;
	}

	int getMatchLastResPos() const {
		return match_last_res_pos_;
	}

	double getPepCTermShift() const {
		return pep_C_term_shift_;
	}

	const PtmPtr& getPepNTermAllowMod() const {
		return pep_N_term_allow_mod_;
	}

	void setPepNTermAllowMod(const PtmPtr& pepNTermAllowMod) {
		pep_N_term_allow_mod_ = pepNTermAllowMod;
	}

	double getPepNTermShift() const {
		return pep_N_term_shift_;
	}

	double getProtCTermShift() const {
		return prot_C_term_shift_;
	}

	double getProtNTermShift() const {
		return prot_N_term_shift_;
	}

	int getTruncLastResPos() const {
		return trunc_last_res_pos_;
	}

	const PtmPtr& getPepCTermAllowMod() const {
		return pep_C_term_allow_mod_;
	}

	const ProtModPtr& getProtCTermAllowMod() const {
		return prot_C_term_allow_mod_;
	}

	const ProtModPtr& getProtNTermAllowMod() const {
		return prot_N_term_allow_mod_;
	}

	bool isAllowPepNMod(){
		return pep_N_term_allow_mod_ == nullptr;
	}

	bool isAllowPepCMod(){
		return pep_C_term_allow_mod_ == nullptr;
	}

	bool isAllowProtNMod(){
		return prot_N_term_allow_mod_ ==nullptr;
	}

	bool isAllowProtCMod(){
		return prot_C_term_allow_mod_ ==nullptr;
	}

	const TruncPtr& getProtCTermAllowTrunc() const {
		return prot_C_term_allow_trunc_;
	}

	const TruncPtr& getProtNTermAllowTrunc() const {
		return prot_N_term_allow_trunc_;
	}

	bool isNTrunc(){
		return n_trunc_;
	}

	bool isCTrunc(){
		return c_trunc_;
	}

	void setPepCTermShift(double pepCTermShift) {
		pep_C_term_shift_ = pepCTermShift;
	}

	void setPepNTermShift(double pepNTermShift) {
		pep_N_term_shift_ = pepNTermShift;
	}

	void setProtCTermShift(double protCTermShift) {
		prot_C_term_shift_ = protCTermShift;
	}

	void setProtNTermShift(double protNTermShift) {
		prot_N_term_shift_ = protNTermShift;
	}

	void setTruncFirstResPos(int truncFirstResPos) {
		trunc_first_res_pos_ = truncFirstResPos;
	}

	void setTruncLastResPos(int truncLastResPos) {
		trunc_last_res_pos_ = truncLastResPos;
	}

	void setMatchFirstResPos(int matchFirstResPos) {
		match_first_res_pos_ = matchFirstResPos;
	}

	void setMatchLastResPos(int matchLastResPos) {
		match_last_res_pos_ = matchLastResPos;
	}

	void setPepCTermAllowMod(const PtmPtr& pepCTermAllowMod) {
		pep_C_term_allow_mod_ = pepCTermAllowMod;
	}

	void setProtCTermAllowMod(const ProtModPtr& protCTermAllowMod) {
		prot_C_term_allow_mod_ = protCTermAllowMod;
	}

	void setProtCTermAllowTrunc(const TruncPtr& protCTermAllowTrunc) {
		prot_C_term_allow_trunc_ = protCTermAllowTrunc;
	}

	void setProtNTermAllowMod(const ProtModPtr& protNTermAllowMod) {
		prot_N_term_allow_mod_ = protNTermAllowMod;
	}

	void setProtNTermAllowTrunc(const TruncPtr& protNTermAllowTrunc) {
		prot_N_term_allow_trunc_ = protNTermAllowTrunc;
	}

	void setNTrunc(bool n_trunc){n_trunc_ = n_trunc;}
	void setCTrunc(bool c_trunc){c_trunc_ = c_trunc;}
	bool isNStrict(){return n_strict_;}
	bool isCStrict(){return c_strict_;}
	void setAlignPrefix(bool is_prefix){is_align_prefix_ = is_prefix;}
	void setAlignsuffix(bool is_suffix){is_align_suffix_ = is_suffix;}
	bool isAlignPrefix(){return is_align_prefix_;}
	bool isAlignSuffix(){return is_align_suffix_;}

	int getId() const {
		return id_;
	}

	void setId(int id) {
		id_ = id;
	}

private:
	int id_=0;
	bool n_trunc_ = false;
	bool n_strict_ = false;
	int trunc_first_res_pos_=0;
	int match_first_res_pos_=0;
	double prot_N_term_shift_;
	double pep_N_term_shift_=0.0;

	TruncPtr prot_N_term_allow_trunc_;
	ProtModPtr prot_N_term_allow_mod_;
	PtmPtr pep_N_term_allow_mod_;

	bool is_align_prefix_ = false;
	bool c_trunc_ = false;
	bool c_strict_ = false;
	int trunc_last_res_pos_=0;
	int match_last_res_pos_=0;
	double prot_C_term_shift_=0.0;
	double pep_C_term_shift_=0.0;

	TruncPtr prot_C_term_allow_trunc_;
	ProtModPtr prot_C_term_allow_mod_;
	PtmPtr pep_C_term_allow_mod_;

	bool is_align_suffix_=false;
};
DiagonalHeaderPtr getShift(DiagonalHeaderPtr shift,int bgn,int end);
DiagonalHeaderPtrVec getNTermShiftListCommon(std::vector<double> best_shifts);
DiagonalHeaderPtrVec getNTermShiftListCompLeft(ProteoformPtr seq,PtmMngPtr mng);
DiagonalHeaderPtrVec getNTermShiftListCompRight(ProteoformPtr seq,PtmMngPtr mng);
void setPrefixSuffix(DiagonalHeaderPtr &header,double c_shift,ProteoformPtr seq,PtmMngPtr mng);
void setProtTermMod(DiagonalHeaderPtr &header,ProteoformPtr seq,PtmMngPtr mng);
void setProtTermTrunc(DiagonalHeaderPtr &header,ProteoformPtr seq,PtmMngPtr mng);
void setPepTermMode(DiagonalHeaderPtr &header,PtmMngPtr mng);
ProtModPtr findProtTermMod(ProtModPtrVec mods,int trunc_len,ResSeqPtr res_seq,double pep_term_shift,double tolerance);
PtmPtr findPepTermMod(PtmPtrVec mods,double shift,double tolerance);
void setAlignPrefSuffic(DiagonalHeaderPtr &header,PtmMngPtr mng);
DiagonalHeaderPtrVec getNTermShiftListTruncPrefix(ProteoformPtr seq);
DiagonalHeaderPtrVec getNTermShiftListTruncsuffix(PrmMsPtr ms,ProteoformPtr seq);
DiagonalHeaderPtrVec get1dHeaders(DiagonalHeaderPtrVec2D headers);
} /* namespace prot */

#endif /* DIAGONAL_HEADER_HPP_ */
