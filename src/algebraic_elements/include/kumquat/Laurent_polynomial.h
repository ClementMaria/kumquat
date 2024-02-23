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

#ifndef KUMQUAT_LAURENT_POLYNOMIAL_H_ 
#define KUMQUAT_LAURENT_POLYNOMIAL_H_

#include <string>
#include <vector>
#include <list>
#include <numeric>
#include <boost/multiprecision/gmp.hpp>
#include <kumquat/number_theory.h>
#include <sstream>

namespace kumquat {

/** \class Laurent_polynomial Laurent_polynomial.h kumquat/Laurent_polynomial.h 
  * \brief Ring of Laurent polynomials with multiprecision integer coefficients.
  *
  * \implements Ring
  * 
  * template IntegerNumber implements IntegerNumber
  * template RingCoefficientStructure implements Ring
  */
template< typename IntegerNumber, 
          typename RingCoefficientStructure >
class Laurent_polynomial {
public:
  static const bool abelian_group = true;
  static const bool pseudo_ring = true;//multiplication undefined a/b*c/d != a/b*(c+d)/d
  static const bool ring = true;
  static const bool principal_ideal_domain = true;
  static const bool field = false;

  typedef RingCoefficientStructure Ring_struct;
  typedef typename Ring_struct::Element Ring_coeff;

/** \brief Initialization of ring.
  **/
  Laurent_polynomial(Ring_struct rstruct) : ring_(rstruct) {}


/** \name Methods for Abelian groups. Implements AbelianGroup.
 * @{ */

/** \brief A (multiprecision) signed integer type, notably for the outer product 
  * encoding the \f$Z\f$-module structure of Abelian groups.
  * 
  * int must be convertible to Integer. Integer must be compatible with the usual 
  * operators +, -, *, / etc. 
  * */
  typedef IntegerNumber Integer;

/** A monomial in the Laurent polynomial, given by a pair (n,c) for encoding a monomial c * q^n.*/
  typedef std::pair<int, Ring_coeff> Monomial;

/** \brief The type of elements of the group.
  * 
  * Elements are represented by a list of non-trivial monomials, ordered by increasing degree.
  * 
  * The additive inverse is the empty list.
  * 
  * Must be copiable. 
  * */
  typedef std::list< Monomial > Element;
/** \brief Return the additive identity 0.
 * 
 * The additive identity can be represented by any pair of integers (0,y) with 
 * \f$y>0\f$.
 * */
  Element additive_identity() { return Element(); }
/** \brief Set \f$a \leftarrow (a+b)\f$. */  
  void plus_equal(Element & a, Element b) { 
    auto ita = a.begin();
    auto itb = b.begin();
    while(ita != a.end() && itb != b.end()) {
      if(ita->first < itb->first ) { ++ita; }
      else {
        if(ita->first > itb->first) {
          a.insert(ita,*itb);
          ++itb;
        }
        else {//ita->first == itb.first
          ring_.plus_equal(ita->second,itb->second);
          auto tmp_ita=ita;
          ++tmp_ita;
          if(ring_.trivial(ita->second)) { a.erase(ita); }
          ita = tmp_ita;
          ++itb;
        }
      }
    }
    while(itb != b.end()) { a.push_back(*itb); ++itb; }
  }
/** \brief Return \f$a+b\f$.*/
  Element plus(Element a, Element b) { 
    plus_equal(a,b);
    return a;
  }

/** \brief Return true iff the element is the additive identity.*/
  bool trivial(Element & a) { return a.empty(); }

/** \brief Set \f$a \leftarrow z \times a\f$, using the \f$Z\f$-module structure of 
 * the group.
 * */
  void times_equal(Element & a, Integer z) { 
    if(z == 0) { a.clear(); }
    else {
      for(auto & monom : a) {
        ring_.times_equal(monom.second,z);
      }
    }
  }
/** \brief Return \f$ z \times a\f$, using the \f$Z\f$-module structure of the 
 * group.
 * */
  Element times(Element a, Integer z) { 
    times_equal(a,z); 
    return a;
  }
/** \brief Return the additive inverse (-a) of an element a. */
  Element additive_inverse(Element a) { 
    if(trivial(a)) return a;
    return times(a,-1);
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
    std::cout << "Should depend on the field of coefficients.\n";
    if(trivial(a)) { return 0; }
    return -1;
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
    auto ita = a.begin();    auto itb = b.begin();
    while(ita != a.end() && itb != b.end()) {
      if( ita->first != itb->first 
         || !ring_.equal(ita->second,itb->second) ) { return false; }
      else { ++ita; ++itb; }
    }
    return (itb == b.end());
  }

/** Construct a Laurent polynomial equal to 0.*/
  Element element() {
    return Element();
  }
/** \brief Converts an Integer into an element of the group.
 * 
 * This input integer is unused as all integers are converted to 0.
 */
  Element element(Integer z) { 
    Element a;
    a.emplace_back(0,z);
    return a; 
  }

/** \brief Construct an Laurent polynomial equal to a single monomial c*q^deg.
 * **/
  Element element(int deg, Ring_coeff c) {
    Element a;
    a.emplace_back(deg,c);
    return a;
  }
/** \brief Construct an Laurent polynomial equal to a single monomial monom.
 * **/
  Element element(Monomial monom) {
    return element(monom.first, monom.second);
  }  
/** \brief Return true iff the input a is equal to the additive identity 0.*/
  bool trivial(Integer z) { return z == 0; }
/* @} */  // end AbelianGroup methods

/** \name Methods for pseudo rings. Implements PseudoRing.
 * @{ */

/** \brief Set \f$a \leftarrow z \times a\f$, using the \f$Z\f$-module structure of 
 * the group.
 * */
  void times_equal(Element & a, Monomial monom) { 
    if(monom.second == 0) { a.clear(); return; }
    for(auto & mon : a) {
      mon.first += monom.first;
      ring_.times_equal(mon.second,monom.second);
    }
  }
/** Set a <- (a*b). */
  void times_equal(Element & a, Element b) {
    Element prod_ab = additive_identity();
    for(auto monom_a : a) {
      Element cp_b(b);
      times_equal(cp_b,monom_a);
      plus_equal(prod_ab,cp_b);
    }
  }
// /** Return a*b.*/
  Element times(Element a, Element b) {
    times_equal(a,b);
    return a;
  }

/** Set a <- a^p. */
    void pow_equal(Element & a, Integer p) {
      Element cp_a(a);
      for(int i=0; i<p;++i) { times_equal(a,cp_a); }
    }
/** Return a <- a^p for a positive integer p. */
    Element pow(Element a, Integer p) {
      pow_equal(a,p);
      return a;
    }
/* @} */  // end pseudo ring methods
/** \name Methods for rings. Implements Ring.
 * @{ */
/** Return the multiplicative identity 1.*/
    Element multiplicative_identity() {
      return element(1);
    }
/** Return the multiplicative inverse of a Element a if it exists, and return the additive identity 0 otherwise.*/
    // Element multiplicative_inverse(Element a) {
    //   return additive_identity();
    // }
/* @} */  // end ring methods


    // Ring_coeff evaluate(Element a, Ring_coeff x) {
    //   std::cout << "to do evaluate\n";
    //   return x;
    // }


  std::string to_string(Element a) {
    if(trivial(a)) { return "0"; }
    std::stringstream ss;
    bool first = true;
    for(auto monom : a) {
      if(first) {
        first = false;
        ss << ring_.to_string(monom.second) << ".q^" << monom.first << " ";
      }
      else {
        ss << "+ " << ring_.to_string(monom.second) << ".q^" << monom.first << " ";
      }
    }
    return ss.str();
  }




private:
  Ring_struct ring_;

};


}  //namespace kumquat

#endif //KUMQUAT_LAURENT_POLYNOMIAL_H_


