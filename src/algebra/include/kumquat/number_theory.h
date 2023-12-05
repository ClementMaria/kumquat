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
#include <random>

namespace kumquat {

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

/** \brief Test whether x is a quadratic residue modulo a prime power p^k (including case x=0), assuming gcd(x,p) = 1 (if x not 0), and 0 <= x < p^k.
 * 
 * p must be a prime number, even or odd.*/
template<typename NumberType>
bool quadratic_residue_subroutine(NumberType x, NumberType p, NumberType k) {
  if(x == 0) { return true; }
  //now x != 0
  if(p == 2) { //x is odd to ensure gcd(x,2)=1
    if(k==1) { return true; }//all numbers are quad res mod 2
    else {
      if(k == 2) { return (x%4)==1; }//x odd quad res mod 4 iff x=1 mod 4 (because all odd squares are 1 mod 4)
      if( k > 2) { return (x%8)==1; }//x odd quad res mod 2^k, k>2, iff x=1 mod 8
    } 
  }
  //else p odd prime
  //An arbitrary number x is a non-trivial residue mod an odd prime p iff its legendre symbol is 1
  //A number relatively prime to an odd prime p is a residue modulo any power of p if and only if it is a residue modulo p
  if( x % p == 0) return true;//case trivial
  return legendre_symbol(x,p) == 1;
}
/** \brief Test whether x is a quadratic residue modulo a prime power p^k (including case x=0), assuming gcd(x,p) = 1 (if x not 0), and 0 <= x < p^k.
 * 
 * p must be a prime number, p_pow_k a power of p, (p even or odd).*/
template<typename NumberType>
bool quadratic_residue(NumberType x, NumberType p_pow_k, NumberType p) {
  NumberType k=0;
  while(p_pow_k != 1) { p_pow_k /= p; ++k; }//compute k such that p_pow_k = p^k
  return quadratic_residue_subroutine(x,p,k);
}
/** \brief Compute a solution y to y^2 = x mod p, for an odd prime p, gcd(x,p) = 1, and 0<x<p. The integer x must be a quadratic residue mod p.
 * 
 * Implements the Tonelli-Shanks algorithm for square root modulo an odd prime.
 * */
template<typename NumberType>
NumberType tonelli_shanks(NumberType x, NumberType p) {
  //find q such that (p-1) = q * 2^s with q odd
  NumberType q = p-1;
  NumberType s = 0;
  while((q % 2) == 0) { q /= 2; ++s; }
  //find a z, 0<z<p, such that z is a quadratic residue mod p
  //expected number of tries = 2;
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> distrib_mod_p((int)1, (int)(p-1));
  NumberType z;
  while(true) {
    z = (NumberType)(distrib_mod_p(gen));
    if(!quadratic_residue(z,p,p)) { break; }
  }
  //z is a quadratic residue mod p, and p = q*2^s, q odd
  NumberType m = s;
  NumberType c = pow(z,q,p);
  NumberType t = pow(x,q,p);
  NumberType power_x = q+1;
  power_x /= 2;
  NumberType r = pow(x,power_x,p);
  while(true) {
    if(t==0) {return 0;}
    if(t==1) {return r;}
    //repeated squaring
    NumberType i=0;
    NumberType pow_t=t;//t^{2^i}
    while(pow_t != 1) { pow_t = (pow_t * pow_t)%p; ++i; }
    NumberType power_c = m-i-1;
    power_c = pow((NumberType)2,power_c);
    NumberType b = pow(c, power_c, p );
    m = i % p;
    c = (b*b) % p;
    t = (t * (b*b)) % p;
    r = (r*b) % p;
  }
}
/** \brief Find a solution y to y^2 = x mod p^k. x must be a quadratic residue mod p^k, which can be checked with quadratic_residue.
 * 
 * p must be a prime number, x must be a quadratic residue modulo p^k, and k >= 1.
 * 
 * A number is a quadratic residue modulo an odd prime power p^k iff it is a quadratic resiude mod p iff its Legendre symbol is 1 (excluding the case 0)
 * */
template<typename NumberType>
NumberType solve_quadratic_residue_subroutine(NumberType x, NumberType p, NumberType k) {
  //we assume 0<= x < p^k and x is a quadratic resiude modulo p^k
  if(x == 0) { return 0; } //0^2 = 0 mod p^k
  if(x == 1) { return 1; }
  if(p == 2) {
    //here k must be > 1 as otherwise cases x = 0,1 has been already caught.
    //also case k==2 is done as x is a quad residue mod 4 iff x=1 mod 4
    //hence k>2 and x = 1 mod 8 
    NumberType res=1; //solution to res^2 = x mod 8
    NumberType pow_two = 8;
    for(int j=4 ; j <= k; ++j) {//after iteration j, res*res = x mod 2^j
      //pow_two == 2^{j-1}
      NumberType lambda = (res*res - x) / pow_two; //
      if( (lambda % 2) != 0 ) { res += (pow_two/2); }
      //else res = res 
      pow_two *= 2;
   }
    return res;
  }
  //now p is odd, x!=0,1
  //solve res^2 = x mod p
  NumberType res_p;
  if( (p%4) == 3) { 
    NumberType power_x = p+1;
    power_x /= 4;
    res_p = pow(x, power_x , p); 
  }//4 | p+1 -> legendre==1
  else {//p%4 == 1, use the slower Tonelli-Shanks algorithm
    res_p = tonelli_shanks(x,p);
  }
  //res_p is a solution to res_p^2 = x mod p
  //lift res_p to a solution to res^2 = x mod p^k via Hensel's lemma
  NumberType pow_p = p;
  for(int j=2; j<=k ; ++j) {
    //turn a solution to res^2 = x mod p^{j-1} to a solution res'^2 = x mod p^j
    NumberType twice_res_p = (2 * res_p) % pow_p; 
    NumberType y = inverse( twice_res_p , pow_p); //such that 2*y*res = 1 mod p^j-1
    pow_p *= p;//p^j
    res_p = (res_p - (res_p*res_p - x) * y) % pow_p;
    if(res_p < 0) { res_p += pow_p; }
  }
  return res_p;
}
/** \brief Find a solution y to y^2 = x mod p^k, where p^k is given as a single number p_pow_k. x must be a quadratic residue mod p^k, which can be checked with quadratic_residue.
 * 
 * p_pow_k must be a power of a prime number p (odd or even), x must be a quadratic residue modulo p^k
 * 
 * */
template<typename NumberType>
NumberType solve_quadratic_residue(NumberType x, NumberType p_pow_k, NumberType p) {
  NumberType k=0;
  while(p_pow_k != 1) { p_pow_k /= p; ++k; }//compute k such that p_pow_k = p^k
  return solve_quadratic_residue_subroutine(x,p,k);
}

} //namespace kumquat

#endif //KUMQUAT_NUMBER_THEORY_H_
