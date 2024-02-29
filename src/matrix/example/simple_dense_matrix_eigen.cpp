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

#include <iostream>
#include <kumquat/Z.h>
#include <kumquat/Dense_matrix.h>
#include <kumquat/Dense_matrix_eigen.h>

using namespace kumquat;

// using boost::multiprecision::mpf_float_1000 float_mp;
// using Eigen::Matrix<float_mp, Dynamic, Dynamic> Matrix_mp;

int main() {
	int n = 10;
	Dense_matrix M(n,n,Z());//nxn uninitialized integer matrix
	for(int i=0; i<n; ++i) {
		for(int j=0; j<n; ++j) {
			M(i,j) = i*j + 1;//symmetric
		}
	}
	//convert into a a Dense_matrix_eigen
	Dense_matrix_eigen Meig(M);
	auto eigenval = Meig.eigenvalues();
	std::sort(eigenval.begin(),eigenval.end());
	for(auto lambda : eigenval) { std::cout << lambda << " "; }
	std::cout << std::endl; 

	return 0;
}