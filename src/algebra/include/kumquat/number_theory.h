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


#ifndef KUMQUAT_NUMBER_THEORY_H_ 
#define KUMQUAT_NUMBER_THEORY_H_

#include <boost/integer/extended_euclidean.hpp>//boost extended gcd
#include <boost/integer/common_factor.hpp>

namespace kumquat {

/** \brief Fast exponentiation of the represenetd integer x, modulo p.
 * 
 * param[in] base An arbitrary number.
 * param[in] exp A positive integer.
 * param[in] mod A non-negative integer. If 0, we compute the power (without mod)
 * @return the integer \f$base^{\operatorname{exp}} \mod \operatorname{mod}\f$.
 * */
template<typename NumberType>
NumberType pow(NumberType base, NumberType exp, NumberType mod = 0) {
  NumberType res = 1;
  if(mod != 0) {
    while (exp > 0) {
      if (exp & 1) { res = (res * base) % mod; }
      base = (base * base) % mod;
      exp >>= 1;
    }
  }
  else {
    while (exp > 0) {
      if (exp & 1) { res = (res * base); }
      base = (base * base);
      exp >>= 1;
    }
  }
  return res;
}
/** \brief Return the Legendre symbol of the represented integer x, w.r.t. a prime 
  * number p, with the explicit formula \f$x^{(p-1) / 2} \mod p\f$, using fast 
  * exponentiation. 
  * 
  * Takes value in \f$\{-1,0,1\}\f$. p must be an odd prime number.
  * */ 
template<typename NumberType>
NumberType legendre_symbol(NumberType x, NumberType p) {
  return kumquat::pow<NumberType>(x%p,(p-1)/2,p);
}

/** \brief Test whether x is a quadratic residue modulo p.
 * 
 * p must be an odd prime number.*/
template<typename NumberType>
bool quadratic_residue(NumberType x, NumberType p) {
  if( x % p == 0) return true;
  return legendre_symbol(x,p) == 1;
}
/** \brief Find a solution y to y^2 = x mod p.
 * 
 * p must be an odd prime number, x must be a quadratic residue modulo p.*/
template<typename NumberType>
NumberType solve_quadratic_residue(NumberType x, NumberType p) {
  NumberType res = 0;
  NumberType x_mod_p = x % p;
  while( ((res*res) % p) != x_mod_p) { ++res; }
  return res;
}
/** \brief Return the inverse of x mod m if it exists, 0 otherwise.
 * 
 * Always returns a non-negative integer y, 0<= y <m*/
template<typename NumberType>
NumberType inverse(NumberType x, NumberType m) {
  auto bezout = boost::integer::extended_euclidean(x, m);
  if(bezout.gcd != 1) { return 0; }
  auto y = bezout.x % m;
  if(y < 0) { return m+y; } 
  return y;
}

/** \brief Compute the division x/y in the PID.
 * 
 * Return the value q such that x = q*y + r, with 0 \leq r < |y|.*/ 
template<typename NumberType>
NumberType division(NumberType x, NumberType y) {
  return x/y;
}
/** \brief Compute the remainder of the division x/y in the PID.
 * 
 * Return the value r such that x = q*y + r, with 0 \leq r < |y|.
 * */ 
template<typename NumberType>
NumberType remainder(NumberType x, NumberType y) {
  return x%y;
}

template<typename NumberType>
NumberType times(NumberType x, NumberType y) {
  return x*y;
}

/** Return the gcd of x and y.*/
template<typename NumberType>
NumberType gcd(NumberType x, NumberType y) {
  return boost::integer::gcd(x, y);
}
/** Return x/gcd(x,y).*/
template<typename NumberType>
NumberType gcd_complement(NumberType x, NumberType y) {
  return (x/gcd(x,y));
}
/** \brief Compute the extended greatest common divisor of two elements of 
 * the ring. 
 * 
 * Return a triple (u,v,gcd) opf ring elements such that gcd is the greatest 
 * common divisor of x and y, and (u,v) satiosfies the Bezout identity:
 * u*x + v*y = gcd, for + the ring addition and * the ring multiplication. 
 * */ 
template<typename NumberType>
std::tuple<NumberType,NumberType,NumberType> extended_gcd(NumberType x, 
                                                          NumberType y) {
  auto res_boost = boost::integer::extended_euclidean(x, y);
  std::tuple<NumberType,NumberType,NumberType> res(res_boost.x, res_boost.y, 
                                                   res_boost.gcd); 
  return res;
}

} //namespace kumquat

#endif //KUMQUAT_NUMBER_THEORY_H_
