/*    This file is part of the KUMQUAT Library -
 *    https://kumquat.inria.fr/ 
 *    - which is a licence protected library. See file LICENSE 
 *    or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef KUMQUAT_DENSE_MATRIX_EIGEN_H_ 
#define KUMQUAT_DENSE_MATRIX_EIGEN_H_

#include <iomanip>
#include <vector>
#include <Eigen/Eigenvalues> 
#include <boost/multiprecision/gmp.hpp>


namespace kumquat {

/** \class Dense_matrix_eigen Dense_matrix_eigen.h kumquat/Dense_matrix_eigen.h 
 * \brief An interface to the dense matrix type of the library Eigen, with real multiprecision coefficients. */
// template< class CoefficientStructure >
class Dense_matrix_eigen : public ScalarRingOperations, public ObjectMonoidalCategory {
public:

	typedef boost::multiprecision::mpf_float_1000 float_mp;
	typedef Eigen::Matrix<float_mp, Dynamic, Dynamic> Matrix_mp;

/** \brief Copy constructor from another matrix type with accessor. Includes the copy constructor.*/ 
	template<typename MatrixType>
	Dense_matrix_eigen(const MatrixType& other) {
		mat_(other.num_rows(), other.num_columns());
		for(int i=0; i<other.num_rows(); ++i) {
			for(int j=0; j<other.num_columns(); ++j) {
				mat_(i,j) = other(i,j);
			}
		}
	}
/** \brief Return the total number of rows in the matrix.*/
	size_t num_rows() { return mat_.rows(); }
/** \brief Return the total number of columns in the matrix.*/
	size_t num_columns() { return mat_.cols(); }
/** \brief Return a vector of all the eigenvalues of the matrix.*/
	std::vector<float_mp> eigenvalues() {
		EigenSolver<Matrix_mp> es(mat_,false);
		std::vector<float_mp> res(es.eigenvalues().begin(), es.eigenvalues().end());
		return res;
	}

private:
	Matrix_mp mat_;

};

} // namespace kumquat