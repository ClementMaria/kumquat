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
  Z PID();
  Dense_matrix< Z > mat(5,7,PID);

  std::cout << "- Fill the matrix with mat[i][j] == i*j+1.\n";
  auto product_plus_one = [&](size_t i, size_t j) -> Z::Element { 
    return PID.element(i * j + 1); 
  };

  mat.fill(product_plus_one);
  std::cout << mat;

  mat.smith_normal_form()

}
