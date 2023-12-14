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
  {//odd p
    size_t n = 10;
    int p = 5;
    int max_pow = 7;
    std::cout << "*** case p = " << p << "\n";
    Q_mod_Z_mp::Integer p_mp = p;
    Q_mod_Z_mp::Integer max_pow_mp = max_pow;
    Q_mod_Z_mp qmodz(p);//group Q_5/Z of fractions a/5^k, 0 <= a < 5^k, k>0
    Dense_matrix< Q_mod_Z_mp > mat(n,n,qmodz);

    std::cout << "- The matrix is " << mat.num_rows() << " by " 
                                    << mat.num_columns() << "\n";

    // std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(314); // mersenne_twister_engine seeded with rd()
    // std::uniform_int_distribution<Q_mod_Z_mp::Integer> distrib_power(1, max_pow);
    // std::uniform_int_distribution<Q_mod_Z_mp::Integer> distrib_nomi(0, 
    //                                                   kumquat::pow(p_mp, max_pow)-1);

    std::uniform_int_distribution<int> distrib_power(1, max_pow);
    std::uniform_int_distribution<int> distrib_nomi(0, 
                                                      std::pow(p, max_pow)-1);


    std::cout << "- Fill the matrix with mat[i][j] = mat[j][i] = random number.\n";
    auto random_qmodz = [&]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j) -> Q_mod_Z_mp::Element { 
      Q_mod_Z_mp::Integer nom = distrib_nomi(gen);
      Q_mod_Z_mp::Integer k = distrib_power(gen);
      Q_mod_Z_mp::Integer den = kumquat::pow(p_mp, k);
      return qmodz.element(nom,den); 
    };

    for(size_t i=0; i<n; ++i) {
      for(size_t j=i; j<n; ++j) {
        auto x = random_qmodz(i,j);
        mat[i][j] = x;
        mat[j][i] = x;
      }
    }

    std::cout << "The matrix is:\n";
    std::cout << mat << "\n";

    mat.diagonalize_gram_matrix_Qp_mod_Z();

    std::cout << "The matrix after normalization is:\n";
    std::cout << mat << "\n";

    Q_U1< Q_mod_Z_mp::Integer > qu1;
    std::cout << "and its Gauss sum is: " << qu1.to_string((mat.gauss_sum_Qp_mod_Z())) << "\n";
  }


  {//even p = 2
    size_t n = 10;
    int p = 2;
    int max_pow = 5;
    std::cout << "*** case p = " << p << "\n";
    Q_mod_Z_mp::Integer p_mp = p;
    Q_mod_Z_mp::Integer max_pow_mp = max_pow;
    Q_mod_Z_mp qmodz(p);//group Q_2/Z of fractions a/2^k, 0 <= a < 2^k, k>0
    Dense_matrix< Q_mod_Z_mp > mat(n,n,qmodz);

    std::cout << "- The matrix is " << mat.num_rows() << " by " 
                                    << mat.num_columns() << "\n";

    // std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(314); // mersenne_twister_engine seeded with rd()
    // std::uniform_int_distribution<Q_mod_Z_mp::Integer> distrib_power(1, max_pow);
    // std::uniform_int_distribution<Q_mod_Z_mp::Integer> distrib_nomi(0, 
    //                                                   kumquat::pow(p_mp, max_pow)-1);

    std::uniform_int_distribution<int> distrib_power(1, max_pow);
    std::uniform_int_distribution<int> distrib_nomi(0, std::pow(p, max_pow)-1);


    std::cout << "- Fill the matrix with mat[i][j] = mat[j][i] = random number.\n";
    auto random_qmodz = [&]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j) -> Q_mod_Z_mp::Element { 
      Q_mod_Z_mp::Integer nom = distrib_nomi(gen);
      Q_mod_Z_mp::Integer k = distrib_power(gen);
      Q_mod_Z_mp::Integer den = kumquat::pow(p_mp, k);
      return qmodz.element(nom,den); 
    };

    for(size_t i=0; i<n; ++i) {
      for(size_t j=i; j<n; ++j) {
        auto x = random_qmodz(i,j);
        mat[i][j] = x;
        mat[j][i] = x;
      }
    }
    
    std::cout << "The matrix is:\n";
    std::cout << mat << "\n";

    mat.diagonalize_gram_matrix_Qp_mod_Z();

    std::cout << "The matrix after normalization is:\n";
    std::cout << mat << "\n";

    Q_U1< Q_mod_Z_mp::Integer > qu1;
    std::cout << "and its Gauss sum is: " << qu1.to_string((mat.gauss_sum_Qp_mod_Z())) << "\n";
  }
  return 0;
}