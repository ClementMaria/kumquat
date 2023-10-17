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
#include <boost/multiprecision/gmp.hpp>
#include <kumquat/number_theory.h>

using namespace kumquat;

typedef boost::multiprecision::mpz_int Int_mp;

int main(int argc, char * argv[]) {
  Int_mp x(5); 
  Int_mp y,z;

  y = 3;
  z = 12;

  Int_mp t = 11;

  std::cout << "Integer 3 = " << y << "\n";
  std::cout << "Integer 3+12 = " << y+z << "\n";
  std::cout << "Integer 5^3 mod 12 = " << pow(x,y,z) << "\n";
  std::cout << "Legendre symbol (5/3) = " << legendre_symbol(x,y) << "\n";
  std::cout << "Legendre symbol (12/5) = " << legendre_symbol(z,x) << "\n";
  std::cout << "Is 12 a quadratic residue mod 5? " << quadratic_residue(z,x) << "\n";
  std::cout << "Legendre symbol (11/5) = " << legendre_symbol(t,x) << "\n";
  std::cout << "Is 11 a quadratic residue mod 5? " << quadratic_residue(t,x) << "\n";
  std::cout << "Solve x^2 = 11 mod 5: x = " << solve_quadratic_residue(t,x) << "\n"; 
}