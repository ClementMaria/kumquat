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
#include <boost/multiprecision/gmp.hpp>
#include <kumquat/number_theory.h>

using namespace kumquat;

typedef boost::multiprecision::mpz_int Int_mp;

int main() {
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
  // std::cout << "Is 12 a quadratic residue mod 5? " << quadratic_residue(z,x) << "\n";
  std::cout << "Legendre symbol (11/5) = " << legendre_symbol(t,x) << "\n";
  // std::cout << "Is 11 a quadratic residue mod 5? " << quadratic_residue(t,x) << "\n";
  // std::cout << "Solve x^2 = 11 mod 5: x = " << solve_quadratic_residue(t,x,5) << "\n"; 


{
  Int_mp prime = 7; //case p%4 == 3
  Int_mp prime_power = 5;
  Int_mp prime_to_pow = pow(prime,prime_power);
  std::cout << "The quadratic residues modulo " << prime_to_pow << ":\n";
  for(int i=0; i<prime_to_pow; ++i) {
    Int_mp i_mp = i;
    if((i_mp % prime) != 0) {
      if(quadratic_residue(i_mp, prime_to_pow, prime)) {
        Int_mp x = solve_quadratic_residue(i_mp, prime_to_pow, prime);
        std::cout << "  " << x << "^2 = " << i << " mod " << prime_to_pow << "\n";
        if( (x*x) % prime_to_pow != i_mp) { std::cout << "Error.\n";}
      }
    }
  }
}

{
  Int_mp prime = 5; //case p%4 == 1
  Int_mp prime_power = 4;
  Int_mp prime_to_pow = pow(prime,prime_power);
  std::cout << "The quadratic residues modulo " << prime_to_pow << ":\n";
  for(int i=0; i<prime_to_pow; ++i) {
    Int_mp i_mp = i;
    if((i_mp % prime) != 0) {
      if(quadratic_residue(i_mp, prime_to_pow, prime)) {
        Int_mp x = solve_quadratic_residue(i_mp, prime_to_pow, prime);
        std::cout << "  " << x << "^2 = " << i << " mod " << prime_to_pow << "\n";
        if( (x*x) % prime_to_pow != i_mp) { std::cout << "Error.\n";}
      }
    }
  }
}

{
  Int_mp prime = 2; //case p == 2
  Int_mp prime_power = 5;
  Int_mp prime_to_pow = pow(prime,prime_power);
  std::cout << "The quadratic residues modulo " << prime_to_pow << ":\n";
  for(int i=0; i<prime_to_pow; ++i) {
    Int_mp i_mp = i;
    if((i_mp % prime) != 0) {
      if(quadratic_residue(i_mp, prime_to_pow, prime)) {
        Int_mp x = solve_quadratic_residue(i_mp, prime_to_pow, prime);
        std::cout << "  " << x << "^2 = " << i << " mod " << prime_to_pow << "\n";
        if( (x*x) % prime_to_pow != i_mp) { std::cout << "Error.\n";}
      }
    }
  }
}

  return 0;
}