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
      base = (base * base) % mod;
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


} //namespace kumquat

#endif //KUMQUAT_NUMBER_THEORY_H_
