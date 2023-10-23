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
#include <kumquat/Z.h>

using namespace kumquat;

int main() {
  Z pid_Z;
  Dense_matrix< Z > mat(5,7,pid_Z);

  std::cout << "- The matrix is " << mat.num_rows() << " by " 
                                  << mat.num_columns() << "\n";

  std::cout << "- Fill the matrix with mat[i][j] == i*j+1.\n";
  auto product_plus_one = [&](size_t i, size_t j) -> Z::Element { 
    return pid_Z.element(i * j + 1); 
  };

  mat.fill(product_plus_one);
  std::cout << mat << "\n\n\n";

  std::vector< Dense_matrix<Z>::Elementary_matrix_operation > ops;
  mat.smith_normal_form(ops);

  std::cout << mat;

  std::cout << "Performing operations: \n";
  for(auto & op : ops) { std::cout << op.to_string() << "\n"; }
}
