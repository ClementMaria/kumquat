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

#ifndef KUMQUAT_Z_MOD_NZ_H_ 
#define KUMQUAT_Z_MOD_NZ_H_

#include <vector>
#include <numeric>
#include <kumquat/number_theory.h>

namespace kumquat {

/** \brief Algebraic structure Z/nZ for an arbirtary n>0.
 *
 * The structure is an abelian group for +, a ring for *, and a field if n is prime.
 * 
 * Elements of the ring Z/nZ are represented by integers (int), between 0 and n-1. 
 * If n is prime, all element of the ring are invertible and Z/nZ is a field.
 *
 * \implements Ring
 */
class Z_mod_nZ {
public:
/** An integer type, in particular for the Z-module structure of the underlying abelian group.*/
  typedef int Integer;
/** The type of elements of the algebraic structure (group/ring/field). */
  typedef int Element;
/** \brief Create the ring Z/nZ. */
  Z_mod_nZ(int N) : N_(N) {//brute force computation of the inverses (if they exist)
    if(N <= 0) { std::cerr << "Instanciation of Z/nZ for n <= 0.\n"; return; }
    inverse_ = std::vector<Element>(N_, 0);
    inverse_[0] = 0;
    inverse_[1] = 1;
    for(int i=2; i<N_; ++i) {
      for(int j=2; j<N_; ++j) {
        if((i*j) % N_ == 1) { inverse_[i] = j; break; } 
      }
    }
  }
/** \brief Convert an integer to an element of the field.*/
  Element element(int z) {
    return (z % N_);
  }
/** \brief Set a <- (a+b). */  
  void plus_equal(Element & a, Element b) {
    //here, 0<= a,b <p by definition
    a = (a+b);
    if(a > N_-1) { a -= N_ ; }//a+b >= p
  }
/** \brief Return the sum a+b of two elements.*/
  Element plus(Element a, Element b) const {
    Element c = a+b;
    if(c > N_-1) { return c-N_; }
    return c;
  }
  /** \brief Set a <- a*b for the \f$Z\f$-module structure of the abelian group. */
  // void times_equal(Element & a, Integer z) {
  //   a = (a * z) % N_;
  // }
/** \brief Set a <- a*b for the ring structure. */
  void times_equal(Element & a, Element b) {
    a = (a * b) % N_;
  }
/** \brief Return the product of an elements and an integer with the \f$Z\f$-module strucutre of the abelian group. */
  // Element times(Element a, Integer z) const {
  //   return (a * z) % N_; 
  // }
/** \brief Return the product of two elements. */
  Element times(Element a, Element b) const {
    return (a * b) % N_; 
  }
/** \brief Return the additive inverse \f$-a\f$ of an element \f$a\f$.*/
  Element additive_inverse(Element a) const {
    return (N_-a);
  }
/** \brief Return the multiplicative inverse \f$a^{-1}\f$ of an element \f$a\f$, if it is invertible. Return 0 otherwise. */
  Element multiplicative_inverse(Element a) const {
    return inverse_[a];
  }
/** \brief Return the additive idendity \f$0 \in \mathbb{Z}/p\mathbb{Z}\f$.*/
  Element additive_identity() const {
    return 0;
  }
/** \brief Return the multiplicative identity \f$1 \in \mathbb{Z}/p\mathbb{Z}\f$.*/
  Element multiplicative_identity() const {
    return 1;
  }
// /** \brief Return a^p for a positive p.*/
  Element pow(Element a, Integer p) {
    return kumquat::pow(a, p, N_);
  }
/** \brief Set a <- a^p for a positive p.*/
  void pow_equal(Element & a, Integer p) {
    a = pow(a,p);
  }

/** Check equality.*/
  bool equal(Element a, Element b) { return a==b; }

/** \brief Return the number of elements N of the abelian group \f$(\mathbb{Z}/N\mathbb{Z},+)\f$. */
  Integer order() const { return N_; }
/** \brief Return the order of the subgroup \f$(\langle a\rangle,+)\f$ of \f$(\mathbb{Z}/N\mathbb{Z},+)\f$.*/
  Integer order(Element a) const { return (N_/std::gcd(a,N_));}
/** Return the rank of the group, i.e., the minimal number of generators. */
  Integer rank() const { return 1; } 
/** \brief Return true iff the input a is equal to the additive identity 0.*/
  bool trivial(Element a) { return a == 0; }

private:
  Integer N_;
  std::vector<Element> inverse_;

};

}  //namespace kumquat


#endif //KUMQUAT_Z_MOD_NZ_H_
