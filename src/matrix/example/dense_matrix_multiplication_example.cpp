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

#include <chrono>
#include <kumquat/Z_mod_nZ.h>
#include <kumquat/Dense_matrix.h>

using namespace kumquat;

// typedef Markov_decision::Proba_float P_float;

int main() {
  Z_mod_nZ G(19);

  int size = 100;
  Dense_matrix< Z_mod_nZ > lmat(size,size,G);

  std::cout << "- The matrix is " << lmat.num_rows() << " by " << lmat.num_columns() << "\n";

  std::cout << "- Fill the matrix with mat[i][j] == i*j+1.\n";
  auto product_plus_one = [&](size_t i, size_t j) -> Z_mod_nZ::Element { 
    return G.element(i * j + 1); 
  };
  lmat.fill(product_plus_one);
  std::cout << lmat;

  Dense_matrix< Z_mod_nZ > rmat(size,size,G);
  rmat.fill(product_plus_one);


std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  lmat.rtimes(rmat);

std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

begin = std::chrono::steady_clock::now();
  
  rmat.ltimes(lmat);

end = std::chrono::steady_clock::now();

std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

  return 0;
}
