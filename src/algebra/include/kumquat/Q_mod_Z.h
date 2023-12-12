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

#ifndef KUMQUAT_Q_MOD_Z_H_ 
#define KUMQUAT_Q_MOD_Z_H_

#include <string>
#include <vector>
#include <numeric>
#include <boost/multiprecision/gmp.hpp>
#include <kumquat/number_theory.h>
#include <sstream>

namespace kumquat {

/** \class Q_mod_Z Q_mod_Z.h kumquat/Q_mod_Z.h 
  * \brief Abelian group \f$\mathbb{Q}/\mathbb{Z}\f$, or its subgroup 
  * \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$ for a prime p.
  *
  * The structure is an Abelian group for +.
  * 
  * Elements are represented by pairs of multiprecision non-negative integers (x,y), 
  * representing fractions \f$x/y\f$, with \f$0 \geq x < y\f$ and gcd(x,y) = 1.
  *
  * \implements AbelianGroup
  * 
  * template IntegerNumber implements IntegerNumber
  */
template<typename IntegerNumber>
class Q_mod_Z {
public:
  static const bool abelian_group = true;
  static const bool pseudo_ring = false;//multiplication undefined a/b*c/d != a/b*(c+d)/d
  static const bool ring = false;
  static const bool principal_ideal_domain = false;
  static const bool field = false;

/** \brief Initialisation of the entire group \f$\mathbb{Q}/\mathbb{Z}\f$.
  * 
  * p_ set to -1 indicates that we construct the entire group.
  * */
  Q_mod_Z() : p_(-1) {}
/** \brief Initialisation of a subgroup \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$ of the 
  * group \f$\mathbb{Q}/\mathbb{Z}\f$.
  * 
  * p_ must be a prime number. Elements are all of the form \f$u / p^k\f$, 
  * \f$ k>0 \f$ and \f$ 0 \leq u < p^k \f$ and gcd(u,p) = 1.
  * */
  Q_mod_Z(int p) : p_(p) {}

/** \name Methods for Abelian groups. Implements AbelianGroup.
 * @{ */

/** \brief A (multiprecision) signed integer type, notably for the outer product 
  * encoding the \f$Z\f$-module structure of Abelain groups.
  * 
  * int must be convertible to Integer. Integer must be compatible with the usual 
  * operators +, -, *, / etc. 
  * */
  typedef IntegerNumber Integer;
/** \brief The type of elements of the group.
  * 
  * Elements are represented by a pair of multiprecision integers \f$(x,y)\f$, 
  * representing the fraction \f$x/y\f$, with \f$0 \leq x < y\f$ and \f$gcd(x,y) = 1\f$. 
  * 
  * Must be copiable. 
  * */
  typedef std::pair< Integer, Integer >   Element;
/** \brief Return the additive identity 0.
 * 
 * The additive identity can be represented by any pair of integers (0,y) with 
 * \f$y>0\f$.
 * */
  Element additive_identity() { return Element(0,1); }
/** \brief Set \f$a \leftarrow (a+b)\f$. */  
  void plus_equal(Element & a, Element b) { 
    a.first = a.first*b.second + a.second*b.first;
    a.second *= b.second;  
    normalize(a);
  }
/** \brief Return \f$a+b\f$.*/
  Element plus(Element a, Element b) { 
    Element a_plus_b(a.first*b.second + a.second*b.first, a.second * b.second);
    normalize(a_plus_b);
    return a_plus_b; 
  }
/** \brief Set \f$a \leftarrow z \times a\f$, using the \f$Z\f$-module structure of 
 * the group.
 * */
  void times_equal(Element & a, Integer z) { 
    if(z<0) {
      z = -1 * z;
      a = additive_inverse(a);//inverse
    }
    a.first *= z; 
    normalize(a); 
  }
/** \brief Return \f$ z \times a\f$, using the \f$Z\f$-module structure of the 
 * group.
 * */
  Element times(Element a, Integer z) { 
    if(z<0) {
      z = -1 * z;
      a = additive_inverse(a);//inverse
    }
    a.first *= z; 
    normalize(a); 
    return a;
  }
/** \brief Return the additive inverse (-a) of an element a. */
  Element additive_inverse(Element a) { 
    if(a.first == 0) return additive_identity();
    return Element(a.second - a.first, a.second); //1-a
  }
/** \brief Return the order of the group. 
 * 
 * It will return the number of element in a finite abelian group, or -1 in case the 
 * group is infinite. 
 * */
  Integer order() { return -1; }
/** \brief Return the order of the Element a. 
 * 
 * It will return the order of the group generated by a if it is finite, or -1 if it 
 * is infinite.
 * */
  Integer order(Element a) { 
    if(a.first == 0) { return 0; }
    auto order_mp = kumquat::gcd_complement(a.second, a.first);
    return order_mp; 
  }
/** \brief Return the rank of the group, i.e., the minimal number of generators. 
 * 
 * Return -1 if the rank is infinite.
 * */
  Integer rank() { return -1; } 
/** Check for equality. Note that the element 0 is ambiguous, as it can be 
 * represented by any fraction /f$0/y/f$ with for any /f$y > 0/f$.
 * */
  bool equal(Element a, Element b) {
    if(a.first == 0) { return b.first == 0; }
    return (a.first * b.second) == (b.first * a.second);
  }
/** \brief Converts an Integer into an element of the group.
 * 
 * This input integer is unused as all integers are converted to 0.
 */
  Element element([[maybe_unused]] Integer z) { return Element(0,1); }
/** \brief Return true iff the input a is equal to the additive identity 0.*/
  bool trivial(Element a) { return a.first == 0; }
  bool trivial(Integer z) { return z == 0; }
/* @} */  // end AbelianGroup methods
/** \name Methods for pseudo rings. Implements PseudoRing.
 * @{ */

/** Set a <- (a*b). */
    // void times_equal(Element & a, Element b) {
    //   a.first = a.first * b.first;
    //   a.second = a.second * b.second;
    //   normalize(a);
    // }
/** Return a*b.*/
    // Element times(Element a, Element b) {
    //   a.first = a.first * b.first;
    //   a.second = a.second * b.second;
    //   normalize(a);
    //   return a;
    // }
/** Set a <- a^p. */
    void pow_equal(Element & a, Integer p) {
      a.first = kumquat::pow(a.first,p);
      a.second = kumquat::pow(a.second,p);
    }
/** Return a <- a^p for a positive integer p. */
    Element pow(Element a, Integer p) {
      a.first = kumquat::pow(a.first,p);
      a.second = kumquat::pow(a.second,p);
      return a;
    }
/* @} */  // end pseudo ring methods


/** \name Methods specific to Q / Z and Q_{p} / Z
 * @{ */
/** \brief Construct the Element \f$x/y \mod \mathbb{Z}\f$ from two integers 
 * /f$x/f$ and /f$y/f$.
 * */ 
  Element element(Integer x, Integer y) { 
    Element z(x,y);
    normalize(z);
    return z; 
  }
/** \brief Return the Integer p in case we represent the subgroup 
 * \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$ of \f$\mathbb{Q}/\mathbb{Z}\f$.
 * 
 * Must return a prime number p, or -1 in case we represent the entire group 
 * \f$\mathbb{Q}/\mathbb{Z}\f$.
 * */
  Integer p() { return p_; }
/** \brief Normnalize a fraction u/v into x/y such that:
 * u/v = x/y modulo Z, and 
 * gcd(x,y) = 1, and 
 * 0 <= x < y.
 */ 
  void normalize(Element & a) {
    normalize_mod_Z(a);
    normalize_fraction(a);
  }
/** \brief Return the numerator of a fraction.*/
  Integer numerator(Element a) { return a.first; }
/** \brief Return the denominator of a fraction.*/
  Integer denominator(Element a) { return a.second; }
  
/* @} */  // end methods specific to Q / Z

// /** \brief Write a fraction a/b, with gcd(a,b)=1, as a product a * 1/b, with a and Integer and 1/b an Element.*/
//   std::pair<Integer,Element> externalize(Element & a) {
//     normalize(a);
//     return std::make_pair(a.first, Element((Integer)1, a.second));
//   }

/** \brief Return a string encoding the fraction.*/
  std::string to_string(Element x) {
    return std::to_string(x.first) + "/" + std::to_string(x.second);
  }


private:
  /* Normalize the element a modulo \f$\mathbb{Z}\f$. 
   * 
   * Note that we do not enforce that \f$\gcd(a.first,a.second) == 1\f$.
   */
  void normalize_mod_Z(Element & a) {
    if(a.first >= a.second) { a.first = (a.first) % (a.second); }
  }
  /* Normalize the fraction such that gcd(numerator,denominator) = 1. 
   */
  void normalize_fraction(Element & a) {
    if(a.first == 0) { a.second = 1; return; }
    auto gcd = kumquat::gcd(a.first,a.second);
    if(gcd != 1) {
      a.first /= gcd;
      a.second /= gcd;
    }
  }
  /** A prime number p in case the structure represents the subgroup 
   * \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$ of \f$\mathbb{Q}/\mathbb{Z}\f$. 
   * 
   * Must be a prime number, or -1.
   * */
  Integer p_;

};

/** \brief Write an Element of Q_mod_Z into a stream, as an unnormalized fraction 
 * x / y.
 * */
// template<typename IntegerNumber>
// std::ostream & operator<<( std::ostream & os, 
//                            typename Q_mod_Z<IntegerNumber>::Element & a) {
//       os << a.first << "/" << a.second;
//   return os;
// }
template<typename IntegerNumber>
std::ostream & operator<<( std::ostream & os, 
                           std::pair<IntegerNumber,IntegerNumber> & a) {
  
  std::stringstream ss;//for alignment purpose
  ss << a.first << "/" << a.second;  
  os << ss.str();
  return os;
}

typedef Q_mod_Z<boost::multiprecision::mpz_int> Q_mod_Z_mp;

}  //namespace kumquat

#endif //KUMQUAT_Q_MOD_Z_H_
