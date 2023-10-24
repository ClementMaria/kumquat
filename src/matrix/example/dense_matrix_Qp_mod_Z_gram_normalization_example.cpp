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
#include <random>
#include <kumquat/Q_mod_Z.h>
#include <kumquat/Dense_matrix.h>
#include <kumquat/number_theory.h>

using namespace kumquat;

int main() {
  size_t n = 7;
  int p = 5;
  Q_mod_Z::Integer p_mp = p;
  Q_mod_Z::Integer max_pow = 3;
  Q_mod_Z qmodz(p);//group Q_5/Z of fractions a/5^k, 0 <= a < 5^k, k>0
  Dense_matrix< Q_mod_Z > mat(n,n,qmodz);

  std::cout << "- The matrix is " << mat.num_rows() << " by " 
                                  << mat.num_columns() << "\n";

  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(314); // mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<Q_mod_Z::Integer> distrib_power(1, max_pow);
  std::uniform_int_distribution<Q_mod_Z::Integer> distrib_nomi(0, 
                                                    kumquat::pow(p_mp, max_pow)-1);

  std::cout << "- Fill the matrix with mat[i][j] = mat[j][i] = random number.\n";
  auto random_qmodz = [&]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j) -> Q_mod_Z::Element { 
    Q_mod_Z::Integer nom = distrib_nomi(gen);
    Q_mod_Z::Integer k = distrib_power(gen);
    Q_mod_Z::Integer den = kumquat::pow(p_mp, k);
    return qmodz.element(nom,den); 
  };

  for(size_t i=0; i<n; ++i) {
    for(size_t j=i; j<n; ++j) {
      auto x = random_qmodz(i,j);
      mat[i][j] = x;
      mat[j][i] = x;
    }
  }

  std::cout << mat << "\n\n\n";

  mat.diagonalize_gram_matrix_Qp_mod_Z_p_odd();

  std::cout << mat << "\n\n\n";

  return 0;
}