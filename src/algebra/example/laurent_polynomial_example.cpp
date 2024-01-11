/*    This file is part of the KUMQUAT Library - https://kumquat.inria.fr/ 
 *    - released under  GNU GENERAL PUBLIC LICENSE v.3
 *    See file LICENSE or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <kumquat/Laurent_polynomial.h>
#include <kumquat/Q.h>

using namespace kumquat;

//Laurent polynomials with rational coefficients and multiprecision integer type
using Integer_t = boost::multiprecision::mpz_int; 
using Ring_coeff = Q< Integer_t >;
using Fraction = Ring_coeff::Element;
using Laurent_poly = Laurent_polynomial< Integer_t, Ring_coeff >;
using L_poly = Laurent_poly::Element;

int main() {
  Ring_coeff Q_coeff;
  Laurent_poly L(Q_coeff);

  Fraction x = Q_coeff.element(-7,2);
  Fraction y = Q_coeff.element(5,31);
  Fraction z = Q_coeff.element(12,3);
  Fraction minus_z = Q_coeff.element(12,-3);
  L_poly a = L.element(-2,z);//monomial 12/3 * q^-2
  L_poly b = L.element(1,x);//monomial -7/2 * q^1
  L_poly c = L.element(1,y);//monomial 5/31 * q^1
  L_poly d = L.element(-2, minus_z);//monomial -12/3 q^-2
  L_poly e = L.element(0,x); //monomial x * q^0 

  std::cout << "a = " << L.to_string(a) << "\n";
  std::cout << "b = " << L.to_string(b) << "\n";
  std::cout << "c = " << L.to_string(c) << "\n";
  L_poly sum = L.plus(b,c);
  std::cout << "b+c = " << L.to_string(sum) << "\n";
  L.plus_equal(sum,a);
  std::cout << "a+b+c = " << L.to_string(sum) << "\n";
  L.plus_equal(sum,d);
  std::cout << "a+b+c+(-a) = " << L.to_string(sum) << "\n";

  Integer_t p = 10;
  std::cout << "10*(b+c) = " << L.to_string(L.times(sum,p)) << "\n";
  L.times_equal(sum,p);
  std::cout << "10*(b+c) = " << L.to_string(sum) << "\n";
  std::cout << "-10*(b+c) = " << L.to_string(L.additive_inverse(sum)) << "\n";

  std::cout << "10*(b+c) == a? " << L.equal(sum,a) << "\n";

  std::cout << "Additive identity = " << L.to_string(L.additive_identity()) << "\n";

  return 0;
}