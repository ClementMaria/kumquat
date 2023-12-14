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

#ifndef KUMQUAT_Q_U1_H_ 
#define KUMQUAT_Q_U1_H_

#include <string>
#include <vector>
#include <numeric>
#include <boost/multiprecision/gmp.hpp>
#include <kumquat/number_theory.h>
#include <sstream>
#include <kumquat/Q.h> 
#include <kumquat/Q_mod_Z.h> 

namespace kumquat {

/** \class Q_U1 Q_U1.h kumquat/Q_U1.h 
  * \brief Field of elements of the form 
  * \f$\frac{a}{b} \times \operatorname{exp}(2 \pi i \frac{r}{s})\f$ for arbitrary 
  * integers \f$a,b,r,s\f$.
  *
  * The structure is an Abelian group for *.
  * 
  * Elements are represented by 4-tuples of multiprecision 
  * integers (a,b,r,s)
  * 
  * \implements AbelianGroup
  * 
  * template IntegerNumber implements IntegerNumber
  */
template<typename IntegerNumber>
class Q_U1 {
public:
  static const bool abelian_group = true;
  static const bool pseudo_ring = false;//multiplication undefined a/b*c/d != a/b*(c+d)/d
  static const bool ring = false;
  static const bool principal_ideal_domain = false;
  static const bool field = false;

  Q_U1() : Q_(), Q_mod_Z_() {};


/** \brief A (multiprecision) signed integer type, notably for the outer product 
  * encoding the \f$Z\f$-module structure of Abelian groups.
  * 
  * int must be convertible to Integer. Integer must be compatible with the usual 
  * operators +, -, *, / etc. 
  * */
  typedef IntegerNumber Integer;  

  typedef typename Q<Integer>::Element       Rational;
  typedef typename Q_mod_Z<Integer>::Element Rational_mod_Z;

/** \brief The type of elements of the group.
  * 
  * Elements are represented by pairs of fractions \f$a/b, r/s\f$ to represent the 
  * number \f$\frac{a}{b} \times \operatorname{exp}(2 \pi i \frac{r}{s})\f$.
  * 
  * The fraction a/b is represented by an element of the rational Q, satisfying 
  * gcd(a,b)=1 and 0 < b.
  * 
  * The fraction r/s is represented by an element of the rational mod the integer 
  * Q_mod_Z, satisfying gcd(r,s)=1 and 0 <= r < s.
  * 
  * In consequence, the elements of Q_U1 are represented by a unique 4-tuple of 
  * integers (a,b,r,s).
  *  
  * Must be copiable. 
  * */
  typedef std::pair< Rational, Rational_mod_Z >   Element;

/** Return the identity for the group law *, as 1/0 exp(2 pi i 0/1 )*/
  Element additive_identity() { 
    return std::make_pair(Q_.multiplicative_identity(), Q_.additive_identity()); 
  }
/** Set a <- (a*b), the group law. */  
  void plus_equal(Element & a, Element b) {
    Q_.times_equal(a.first,b.first);
    Q_mod_Z_.plus_equal(a.second,b.second);
  }
/** Return (a*b), the group law. */  
  void plus(Element a, Element b) {
    Q_.times_equal(a.first,b.first);
    Q_mod_Z_.plus_equal(a.second,b.second);
    return a;
  }
/** Set a<- z*a using the Z-module structure of the group.*/
  void times_equal(Element a, Integer z) {
    Q_.times_equal(a.first,z);
  }
/** Return z*a using the Z-module structure of the group.*/
  Element times(Element a, Integer z) {
    Q_.times_equal(a.first,z);
    return a;    
  }
/** Return the additive inverse (-a) of element a. */
  Element additive_inverse(Element a) {
   return std::make_pair( Q_.multiplicative_inverse(a.first), 
                          Q_mod_Z_.additive_inverse(a.second) );
  }
/** Return the order of the group. It will return the number of element in a finite abelian group, or -1 in case the group is infinite. */
  Integer order() { return -1; }
/** Return the order of the Element a. It will return the order of the group generated by a if it is finite, or -1 if it is infinite.*/
  Integer order(Element a) { 
    if(Q_.trivial(a.first)) { return 0; }    
    if(Q_.equal(a.first, Q_.multiplicative_identity())) {
      return Q_mod_Z_.denominator(a.second);
    }
    return -1;
  }
/** Return the rank of the group, i.e., the minimal number of generators. */
  Integer rank() { return -1; } 
/** Return true iff a and b represent the same element of the group.*/
  bool equal(Element a, Element b) {
    return (Q_.equal(a.first,b.first)) && (Q_mod_Z_.equal(a.second,b.second));
  }
/** Convert an integer to an element of the group, equal to z*1.*/
  Element element(Integer z) {
    return std::make_pair(Q_.element(z,1), Q_mod_Z_.additive_identity());
  }
/** Construct  a/b*e(r/s).*/
  Element element(Integer a, Integer b, Integer r, Integer s) {
    return std::make_pair(Q_.element(a,b), Q_mod_Z_.element(r,s));
  }
/** Check whether an element is trivial.*/
  bool trivial(Element a) {
    return equal(a,additive_identity());
  }
/** Check whether an integer is trivial.*/
  bool trivial(Integer a) { return a == 0; } 
/* @} */  // end AbelianGroup methods

/** \brief Return a string encoding the element.
 * 
 * An element (x,y,r,s) representing the number (x/y)* exp(i 2 pi * (r/s) ) gives the string: x/ye(r/s)
 * */
  std::string to_string(Element x) {
    return Q_.to_string(x.first) + "e(" + Q_mod_Z_.to_string(x.second) + ")";
  }

private:
  Q<Integer>       Q_;
  Q_mod_Z<Integer> Q_mod_Z_;

};

typedef Q_U1<boost::multiprecision::mpz_int> Q_U1_mp;

}  //namespace kumquat

#endif //KUMQUAT_Q_MOD_Z_H_
