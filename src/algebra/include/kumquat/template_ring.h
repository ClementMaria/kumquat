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

#ifndef KUMQUAT_XXX_H_ 
#define KUMQUAT_XXX_H_

#include <vector>
#include <numeric>
#include <kumquat/number_theory.h>
#include <boost/multiprecision/gmp.hpp>//boost multiprecision wrap over gmp
// #include <boost/integer/extended_euclidean.hpp>//boost extended gcd

namespace kumquat {

/** \brief Algebraic structure .
 *
 * The structure is an abelian group for +, and a ring for *.
 * 
 * Elements of the ... are represented by ....
 *
 * \implements Ring
 */
template<typename IntegerNumber>
class X_xx {
public:
  static const bool abelian_group = true;
  static const bool pseudo_ring = true;
  static const bool ring = true;
  static const bool principal_ideal_domain = true;
  static const bool field = false;

  X_xx() {}
/** \name Methods for Abelian groups. Implements AbelianGroup.
 * @{ */
/** \brief An integer type, in particular for the Z-module structure.*/
  typedef IntegerNumber Integer;
/** \brief The type of elements of the ring. Must be copiable. */
  typedef ... Element;
/** \brief Return the additive identity 0.*/
  Element additive_identity() { return Element(0); }
/** \brief Set a <- (a+b). */  
  void plus_equal(Element & a, Element b) { a += b; }
/** Return a+b.*/
  Element plus(Element a, Element b) { return a+b; }
/** \brief Set a<- z*a using the Z-module structure of the group.
 * 
 * Because the types Integer and Element coincide, this is also the product coming from the ring structure of the integers.
 * */
  void times_equal(Element a, Integer z) { a *= Element(z); }
/** \brief Return z*a using the Z-module structure of the group.
 */
  Element times(Element a, Integer z) { return a*Element(z); }
/** Return the additive inverse (-a) of element a. */
   Element additive_inverse(Element a) { return a*Element(-1); }
/** \brief Return the order of the group. 
 * 
 * It will return the number of element in a finite abelian group, or -1 in case the group is infinite. 
 * */
  Integer order() { return ...; }
/** \brief Return the order of the Element a. 
 * 
 * It will return the order of the group generated by a if it is finite, or -1 if it is infinite. */
  Integer order(Element a) { 
    if(trivial(a)) { return 0; }
    return ...; 
  }
/** \brief Return the rank of the group, i.e., the minimal number of generators. */
  Integer rank() { return ...; } 
/** \brief Return true iff a and b represent the same element of the group.*/
  bool equal(Element a, Element b) { return a==b; }
/** Convert an integer to an element of the ring Z.*/
  Element element(Integer z) { return Element(z); }
/** Check whether an element is trivial.*/
  bool trivial(Element a) { return a == Element(0); }

/* @} */  // end AbelianGroup methods



/** \name Methods for rings. Implements Ring.
 * @{ */
/** \brief Return the multiplicative identity 1.*/
  Element multiplicative_identity() { return Element(1); }
/** \brief Set a <- (a*b). */
  void times_equal(Element & a, Element b) { a *= b; }
/** Return a*b.*/
  Element times(Element a, Element b) { return a*b; }
/** Set a <- a^p. */
  void pow_equal(Element & a, Integer p) { a = kumquat::pow(a,p); }
/** Return a <- a^p for a positive integer p. */
  Element pow(Element a, Integer p) { 
    return kumquat::pow(a,p); 
  }
/** \brief Return the multiplicative inverse of an integer a if it exists, and return the additive identity 0 otherwise.*/
  Element multiplicative_inverse(Element a) {
    ...
    return 0;
  }

/* @} */  // end Ring methods

/** \name Methods for principal ideal domains. Implements PrincipalIdealDomain.
 * @{ */

/** \brief Return the absolute value of an integer. 
 * 
 * Plays the role of the Euclidean function for the structure of the PID. See functions division and remainder*/
  Element abs(Element a) { 
    if(a < 0) { times_equal(a,-1); }
    return a;
  }
/** \brief Compute the extended greatest common divisor of two elements of 
 * the ring. 
 * 
 * Return a triple (u,v,gcd) opf ring elements such that gcd is the greatest 
 * common divisor of x and y, and (u,v) satiosfies the Bezout identity:
 * u*x + v*y = gcd, for + the ring addition and * the ring multiplication. 
 * */ 
  std::tuple<Element,Element,Element> extended_gcd(Element x, Element y) {
    return kumquat::extended_gcd(x,y);
  }
/** \brief Compute the division x/y in the PID.
 * 
 * Return the value q such that x = q*y + r, with 0 \leq r < |y|.*/ 
  Element division(Element x, Element y) {
    return x/y;//defined in number_theory.h
  }
/** \brief Compute the remainder of the division x/y in the PID.
 * 
 * Return the value r such that x = q*y + r, with 0 \leq r < |y|.
 * */ 
  Element remainder(Element x, Element y) {
    return x%y;//defined in number_theory.h
  }

/* @} */  // end Principal ideal domain methods


/** \brief Return a string encoding the element.*/
  std::string to_string(Element x) {
    return ...;
  }  

};

}  //namespace kumquat

#endif //KUMQUAT_Z_H_