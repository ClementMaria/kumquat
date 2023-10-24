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

#include <vector>
#include <numeric>
#include <boost/multiprecision/gmp.hpp>

namespace kumquat {

/** \class Q_mod_Z Q_mod_Z.h kumquat/Q_mod_Z.h 
  * \brief Abelian group \f$\mathbb{Q}/\mathbb{Z}\f$, or its subgroup 
  * \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$ for a prime p.
  *
  * The structure is an abelian group for +.
  * 
  * Elements are represented by pairs of multiprecision non-negative integers (x,y), 
  * repesenting fractions \f$x/y\f$, with \f$0 \geq x < y\f$.
  *
  * \implements AbelianGroup
  */
class Q_mod_Z {
public:
/** \brief Initialisation of the entire group \f$\mathbb{Q}/\mathbb{Z}\f$.
  * 
  * p_ set to -1 indicates that we construct the entire group.
  * */
  Q_mod_Z() : p_(-1) {}
/** \brief Initialisation of a subgroup \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$ of the 
  * group \f$\mathbb{Q}/\mathbb{Z}\f$.
  * 
  * p_ must be a prime number. Elements are all of the form \f$u / p^k\f$, 
  * \f$ k>0 \f$ and \f$ 0 \leq u < p^k \f$, but may be represented by an 
  * unnormalized fraction \f$x/y\f$ with \f$\gcd(x,y) \neq 1\f$.
  * */
  Q_mod_Z(int p) : p_(p) {}
/** \brief A (multiprecision) signed integer type, notably for the outer product 
  * encoding the \f$Z\f$-module structure of Abelain groups.
  * 
  * int must be convertible to Integer. Integer must be compatible with the usual 
  * operators +, -, *, / etc. 
  * */
  typedef boost::multiprecision::mpz_int Integer;
/** \brief The type of elements of the group.
  * 
  * Elements are epresented by a pair of multiprecision integers \f$(x,y)\f$, 
  * representing the fraction \f$x/y\f$, with \f$0 \leq x < y\f$. We do not 
  * however enforce that \f$gcd(x,y) = 1\f$. 
  * 
  * Must be copiable. 
  * */
  typedef std::pair< boost::multiprecision::mpz_int,
                     boost::multiprecision::mpz_int >   Element;
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
    normalize_element(a);
  }
/** \brief Return \f$a+b\f$.*/
  Element plus(Element a, Element b) { 
    Element a_plus_b(a.first*b.second + a.second*b.first, a.second * b.second);
    normalize_element(a_plus_b);
    return a_plus_b; 
  }
/** \brief Set \f$a \leftarrow z \times a\f$, using the \f$Z\f$-module structure of 
 * the group.
 * */
  void times_equal(Element & a, Integer z) { 
    a.first *= z; 
    normalize_element(a); 
  }
/** \brief Return \f$ z \times a\f$, using the \f$Z\f$-module structure of the 
 * group.
 * */
  Element times(Element a, Integer z) { 
    a.first *= z; 
    normalize_element(a); 
    return a;
  }
/** \brief Return the additive inverse (-a) of an element a. */
  Element additive_inverse(Element a) { 
    if(a.first == 0) return additive_identity();
    return Element(a.second - a.first, a.second); //1-a
  }
/** Check for equality. Note that the element 0 is ambiguous, as it can be 
 * represented by any fraction /f$0/y/f$ with for any /f$y > 0/f$.
 * */
  bool equal(Element a, Element b) {
    if(a.first == 0) { return b.first == 0; }
    return (a.first * b.second) == (b.first * a.second);
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
    auto order_mp = (a.second / boost::multiprecision::gcd(a.first,a.second));
    return (Integer)order_mp; 
  }
/** \brief Return the rank of the group, i.e., the minimal number of generators. 
 * 
 * Return -1 if the rank is infinite.
 * */
  Integer rank() { return -1; } 

/** \brief Converts an Integer into an element of the group.
 * 
 * This input integer is unused as all integers are converted to 0.
 */
  Element element([[maybe_unused]] Integer z) { return Element(0,1); }

/** \brief Construct the Element \f$x/y \mod \mathbb{Z}\f$ from two integers 
 * /f$x/f$ and /f$y/f$.
 * */ 
  Element element(Integer x, Integer y) { return Element(x % y,y); }

/** \brief Return the Integer p in case we represent the subgroup 
 * \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$ of \f$\mathbb{Q}/\mathbb{Z}\f$.
 * 
 * Must return a prime number p, or -1 in case we represent the entire group 
 * \f$\mathbb{Q}/\mathbb{Z}\f$.
 * */
  Integer p() { return p_; }
/** \brief Return true iff the input a is equal to the additive identity 0.*/
  bool trivial(Element a) { return a.first == 0; }

  void p_normalize(Element & a) {
    normalize_element(a);
    Integer gcd = boost::integer::gcd(a.first,a.second);
    a.first /= gcd;
    a.second /= gcd;
  }

private:
  /** Normalize the element a modulo \f$\mathbb{Z}\f$. 
   * 
   * Note that we do not enforce that \f$\gcd(a.first,a.second) == 1\f$.
   */
  void normalize_element(Element & a) {
    if(a.first >= a.second) { a.first = (a.first) % (a.second); }
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
std::ostream & operator<<(std::ostream & os, Q_mod_Z::Element & a) {
      os << a.first << "/" << a.second;
  return os;
}

}  //namespace kumquat

#endif //KUMQUAT_Q_MOD_Z_H_
