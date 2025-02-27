//
// $Id: MatrixInverse.hpp 4146 2012-11-26 23:46:34Z pcbrefugee $
//
//
// NB: Variations of this file appear in many open source projects,
// with no copyright claims or license statements made in any of them.  
// Assumed to be public domain, or at least as open as the boost license,
// since the farthest back we can seem to trace it is the Boost Wiki at
// http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?Effective_UBLAS/Matrix_Inversion
//

#ifndef _MS_UTIL_MATRIXINVERSE_HPP_
#define _MS_UTIL_MATRIXINVERSE_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace toppic {

namespace matrix_inverse {

// Used in savitzky_golay.cpp

// Matrix inversion routine.
// Uses lu_factorize and lu_substitute in uBLAS to invert a matrix 
bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, 
                  boost::numeric::ublas::matrix<double>& inverse)
{

    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<double> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if( res != 0 ) return false;

    // create identity matrix of "inverse"
    inverse = identity_matrix<double>(A.size1());

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}


/**
* Invert a matrix via gauss-jordan algorithm (PARTIAL PIVOT)
*
* @param m The matrix to invert. Must be square.
* @param singular If the matrix was found to be singular, then this
*        is set to true, else set to false.
* @return If singular is false, then the inverted matrix is returned.
*         Otherwise it contains random values.
*/

//#define T double /// for debug
boost::numeric::ublas::matrix<double>
gjinverse(const boost::numeric::ublas::matrix<double> &m, 
          bool &singular)
{
    using namespace boost::numeric::ublas;

    const size_t size = m.size1();

    // Cannot invert if non-square matrix or 0x0 matrix.
    // Report it as singular in these cases, and return 
    // a 0x0 matrix.
    if (size != m.size2() || size == 0)
    {
        singular = true;
        matrix<double> A(0,0);
        return A;
    }

    // Handle 1x1 matrix edge case as general purpose 
    // inverter below requires 2x2 to function properly.
    if (size == 1)
    {
        matrix<double> A(1, 1);
        if (m(0,0) == 0.0)
        {
            singular = true;
            return A;
        }
        singular = false;
        A(0,0) = 1/m(0,0);
        return A;
    }

    // Create an augmented matrix A to invert. Assign the
    // matrix to be inverted to the left hand side and an
    // identity matrix to the right hand side.
    matrix<double> A(size, 2*size);
    matrix_range<matrix<double> > Aleft(A, 
        range(0, size), 
        range(0, size));
    Aleft = m;
    matrix_range<matrix<double> > Aright(A, 
        range(0, size), 
        range(size, 2*size));
    Aright = identity_matrix<double>(size);

    // Doing partial pivot
    for (size_t k = 0; k < size; k++)
    {
        // Swap rows to eliminate zero diagonal elements.
        for (size_t kk = 0; kk < size; kk++)
        {
            if ( A(kk,kk) == 0 ) // XXX: test for "small" instead
            {
                // Find a row(l) to swap with row(k)
                int l = -1;
                for (size_t i = kk+1; i < size; i++) 
                {
                    if ( A(i,kk) != 0 )
                    {
                        l = i; 
                        break;
                    }
                }

                // Swap the rows if found
                if ( l < 0 ) 
                {
                    std::cerr << "Error:" <<  __FUNCTION__ << ":"
                        << "Input matrix is singular, because cannot find"
                        << " a row to swap while eliminating zero-diagonal.";
                    singular = true;
                    return Aleft;
                }
                else 
                {
                    matrix_row<matrix<double> > rowk(A, kk);
                    matrix_row<matrix<double> > rowl(A, l);
                    rowk.swap(rowl);

/*#if defined(DEBUG) || !defined(NDEBUG)
                    std::cerr << __FUNCTION__ << ":"
                        << "Swapped row " << kk << " with row " << l 
                        << ":" << A << "\n";
#endif*/
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////// 
        // normalize the current row
        for (size_t j = k+1; j < 2*size; j++)
            A(k,j) /= A(k,k);
        A(k,k) = 1;

        // normalize other rows
        for (size_t i = 0; i < size; i++)
        {
            if ( i != k )  // other rows  // FIX: PROBLEM HERE
            {
                if ( A(i,k) != 0 )
                {
                    for (size_t j = k+1; j < 2*size; j++)
                        A(i,j) -= A(k,j) * A(i,k);
                    A(i,k) = 0;
                }
            }
        }

/*#if defined(DEBUG) || !defined(NDEBUG)
        std::cerr << __FUNCTION__ << ":"
            << "GJ row " << k << " : " << A << "\n";
#endif*/
    }

    singular = false;
    return Aright;
}

}

}

#endif // _MATRIXINVERSE_HPP_
