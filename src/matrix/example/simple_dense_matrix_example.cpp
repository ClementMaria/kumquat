/*    This file is part of the KUMQUAT Library - https://kumquat.inria.fr/ 
 *    - released under  GNU GENERAL PUBLIC LICENSE v.3
 *    See file LICENSE or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <kumquat/Dense_matrix.h>
#include <kumquat/Z_mod_nZ.h>

using namespace kumquat;

int main() {
	Z_mod_nZ G(19);
	Dense_matrix< Z_mod_nZ > mat(5,7,G);

	std::cout << "- The matrix is " << mat.num_rows() << " by " << mat.num_columns() << "\n";

	std::cout << "- Fill the matrix with mat[i][j] == i*j+1.\n";
	auto product_plus_one = [&](size_t i, size_t j) -> Z_mod_nZ::Element { 
		return G.element(i * j + 1); 
	};

	mat.fill(product_plus_one);
	std::cout << mat;

	std::cout << "- Exchange column 2 with column 6.\n";
	mat.exchange_col(2,6);
	std::cout << mat;

	std::cout << "- Exchange row 0 with row 1.\n";
	mat.exchange_row(0,1);
	std::cout << mat;

	std::cout << "- Set [row 0] <- [row 0] + 18 * [row 1].\n";
	mat.plus_equal_row(0, 1, G.element(18));
	std::cout << mat;

	std::cout << "- Set [col 0] <- [col 0] + 2 * [col 1].\n";
	mat.plus_equal_column(0, 1, G.element(2));
	std::cout << mat;

	std::vector< std::tuple<size_t,size_t,Z_mod_nZ::Element> > col_ops;
	mat.column_echelon_form(col_ops);
	std::cout << mat;

	std::cout << "dimension of kernel = " << mat.dim_kernel() << "  and dimension of image = " << mat.dim_image() << "\n";

}