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

#ifndef KUMQUAT_JONES_POLYNOMIAL_H_ 
#define KUMQUAT_JONES_POLYNOMIAL_H_

#include <kumquat/Rational_function_integral_mp>

namespace kumquat {

/** \brief Class to produce and store the matrices, deduced from appropriate fusion categories, used to defined quantum invariants under several constructions (Reshetikhin-Turaev, Turaev-Viro-Barrett-Westbury, etc).
 *
 * The approach is lazy and keep in memory the useful matrices that have already 
 * been computed.
 **/
// template<int N>
class Quantum_data {
public:
  typedef boost::multiprecision::mpz_int Integer;
  typedef Rational_function_integral_mp Rational_f;

  Quantum_data();

/** Jones polynomials for braid. We follow "Ohtsuki - Quantum invariants: A study of knots, 3-manifolds, and their sets" (book p.77).
 * 
 * Define X = q^{1/2} everywhere.
 * */

/** Types of quantum groups:
 * Uqsl2_genq(C) for \f$U_q(\operatorname{sl}_2(\mathbb{C}))\f$ at a generic q.
 */
enum Quantum_group_type {Uqsl2C_genq};

/** \brief Return the quantum integer 
 * \f$[n]_q = \frac{q^{n/2} - q^{-n/2}}{q^{1/2}-q^{-1/2}}\f = 
 *          = \frac{X^{n}-X^{-n}}{X - X^{-1}}
 *          = \frac{X^{2n}-1}{X^{n+1}-X^{n-1}}\f$
 * 
 * (is also a polynomial \f$q^{1-n}(1+q+q^2+\ldots+q^{n-1})=q^{-n+1}+ ... + q^{-1} + 1\f$).
 */
  Rational_f quantum_integer(int n) {
    if(n<0) { return Rational_f(); }//by convention, [n]=0, n negative
    if(n<2) { return Rational_f(1); }//[0]=[1]=1
    std::vector<Integer> num(2*n+1,0);
    num[0] = -1; num[2*n] = 1;
    std::vector<Integer> den(n+1,0);
    den[n-1] = -1; den[n+1] = 1;
    Rational_f q_int(num,den);//==0
    return q_int;
  }
/** \brief Return the quantum factorial 
 * \f$[n]_q! = [1]_q \cdot [2]_q \ldots [n]_q\f$.
 */
  Laurent_poly quantum_factorial(int n) {
    if(n<0) { return Rational_f(0); }//by convention, [n]!=0, n negative
    if(n<2) { return Rational_f(1); }//[0]!=[1]!=1
    auto q_fact = Rational_f(1);//1 = [1]_q
    for(int i=2; i<n+1; ++i) {
      q_fact *= quantum_integer(i);//*= [i]_q
    }
    return q_fact;
  }

/** \brief Return the q-exponential map of a matrix M, equal to:
 * \f[\operatorname{exp}_q(M) = \displaystyle\sum_{n=0}^{+\infty} \frac{q^{n(n-1)/4}}{[n]_q!} M^n\f].
 * */
  Matrix q_exponential_map(Matrix &M) {
    Matrix powM = 
  }

/** \name Methods for the implementation of framed oriented tangle/link invariants derived from an irreducible representation \f$\rho: A \to \operatorname{End}(V)\f$ of a ribbon Hopf algebras \f$A\f$ (with data) into a vector space \f$V\f$, as described in "Ohtsuki - Quantum invariants: A study of knots, 3-manifolds, and their sets" Chapter 4, using the same notations.
 * @{ */
  Matrix rho_u(int N, Quantum_group_type Uqg) {
    switch(Uqg) {
    case Uqsl2C:
      break;
    default: std::cerr << "Invalid quantum group.\n";
    }
  }
  Matrix rho_v(int N, Quantum_group_type Uqg) {
    switch(Uqg) {
    case Uqsl2C:
      break;
    default: std::cerr << "Invalid quantum group.\n";
    }
  }
  Matrix rho_n(int N, Quantum_group_type Uqg) {
    switch(Uqg) {
    case Uqsl2C:
      break;
    default: std::cerr << "Invalid quantum group.\n";
    }
  }
  Matrix rho_np(int N, Quantum_group_type Uqg) {
    switch(Uqg) {
    case Uqsl2C:
      break;
    default: std::cerr << "Invalid quantum group.\n";
    }
  }

/** \brief Return the R-matrix as defined in Ohtsuki (4.30) p.77.*/
  Matrix h_matrix(int N, Quantum_group_type Uqg) {
    switch(Uqg) {
    case Uqsl2C:
      break;
    default: std::cerr << "Invalid quantum group.\n";
    }
  }
/** \brief Return the R-matrix as defined in Ohtsuki (4.30) p.77.*/
  Matrix R_matrix(int N) {

  }
/** \brief Return the inverse of the R-matrix for the N-th colored Jones polynomial of a braid.*/
  Matrix R_inv_matrix(int N) {

  }


/* @} */  // end AbelianGroup methods


  Laurent_poly jones_polynomial(Braid b, int N = 2) {

  }

private:
  Laurent_poly_struct L_;
  std::map<int, >
};

} //namespace kumquat

#endif // KUMQUAT_JONES_POLYNOMIAL_H_
