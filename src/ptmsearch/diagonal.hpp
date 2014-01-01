/*
 * diagonal.hpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#ifndef DIAGONAL_HPP_
#define DIAGONAL_HPP_

namespace prot {

template <class T>
class Diagonal {
public:
	Diagonal(){};
	Diagonal(DiagonalHeaderPtr header){
		header_ = header;
	};
	Diagonal(DiagonalHeaderPtr header,std::vector<T> pair_ptr_list){
		header_ = header;
		pair_ptr_list_ = pair_ptr_list;
	};

	unsigned int size(){
		return pair_ptr_list_.size();
	}
	DiagonalHeaderPtr getHeader(){
		return header_;
	}
	std::vector<T> getDiagPair(){
		return pair_ptr_list_;
	}

	T getDiagPair(int i){
		return pair_ptr_list_[i];
	}

private:
	DiagonalHeaderPtr header_;
	std::vector<T> pair_ptr_list_;
};

} /* namespace prot */

#endif /* DIAGONAL_HPP_ */
