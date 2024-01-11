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

#ifndef KUMQUAT_JONES_POLYNOMIAL_H_ 
#define KUMQUAT_JONES_POLYNOMIAL_H_

#include <boost/multiprecision/gmp.hpp>

namespace kumquat {

/** Class to handle R-matrices and compute the colored Jones polynomials of links and braids.
 * 
 * 
 **/
class Jones_polynomial {
public:
  typedef boost::multiprecision::mpz_int Integer;
  typedef Laurent_polynomial<Integer, Z<Integer> > Laurent_poly;

  Jones_polynomial();

/** Jones polynomials for braid. We follow "Ohtsuki - Quantum invariants: A study of knots, 3-manifolds, and their sets" (book p.77).*/

/** \brief Return the matrix h for the N-th colored Jones polynomial of a braid.*/
  Matrix h_matrix(int N) {

  }
/** \brief Return the R-matrix for the N-th colored Jones polynomial of a braid.*/
  Matrix R_matrix(int N) {

  }
/** \brief Return the inverse of the R-matrix for the N-th colored Jones polynomial of a braid.*/
  Matrix R_inv_matrix(int N) {

  }



  Laurent_poly jones_polynomial(Braid b, int N = 2) {

  }

private:
};

} //namespace kumquat

#endif // KUMQUAT_JONES_POLYNOMIAL_H_
