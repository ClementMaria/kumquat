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

/** \brief Class implementing a set of quantum data and methods, such as quantum integers, quantum exponential, etc.
 * 
 * 
 * 
 * 
 **/
// template<int N>
class Quantum_data {
public:
  typedef boost::multiprecision::mpz_int Integer;
  typedef Rational_function_integral_mp Rational_f;
  Quantum_data();

/** Jones polynomials for braid. We follow "Ohtsuki - Quantum invariants: A study of knots, 3-manifolds, and their sets" (book p.77).
 * 
 * Define X = q^{1/4} everywhere.
 * */

/** \brief Return \f$q^{1/2}-q^{-1/2} = X^2 - X^{-2} = (X^4-1)/(X^2)\f$ as a rational function.*/
  Rational_f quantum_half() {
    std::vector<Integer> num(5,0); num[0]=-1; num[4]=1;
    std::vector<Integer> den(3,0); den[2]=1;
    return Rational_f(num,den);
  }
/** \brief Return \f$X^k = q^{k/4}\f$ for any integer k.*/
  Rational_f quantum_monomial(int k) {
    if(k>=0) {//return X^k
      std::vector<Integer> num(k+1,0); num[k]=1;
      std::vector<Integer> den(1,1);//==1
      return Rational_f(num,den);
    }
    if(k<0) {
      std::vector<Integer> num(1,1);//==1
      k *= -1;
      std::vector<Integer> den(k+1,0); den[k]=1;
      return Rational_f(num,den);
    }
  }
/** \brief Return the quantum integer 
 * \f$[n]_q = \frac{q^{n/2} - q^{-n/2}}{q^{1/2}-q^{-1/2}}\f = 
 *          = \frac{X^{2n}-X^{-2n}}{X^2 - X^{-2}}
 *          = \frac{X^{4n}-1}{X^{2(n+1)}-X^{2(n-1)}}\f$
 */
  Rational_f quantum_integer(int n) {
    if(n<0) { return Rational_f(); }//by convention, [n]=0, n negative
    if(n<2) { return Rational_f(1); }//[0]=[1]=1
    std::vector<Integer> num(4*n+1,0);
    num[0] = -1; num[4*n] = 1;//X^{4n} - X^0
    std::vector<Integer> den(2*n+3,0);
    den[2*n+2] = 1; den[2*n-2] = -1;
    return Rational_f(num,den);//==0
  }
/** \brief Return the quantum factorial 
 * \f$[n]_q! = [1]_q \cdot [2]_q \ldots [n]_q\f$.
 */
  Rational_f quantum_factorial(int n) {
    if(n<0) { return Rational_f(0); }//by convention, [n]!=0, n negative
    if(n<2) { return Rational_f(1); }//[0]!=[1]!=1
    auto q_fact = Rational_f(1);//1 = [1]_q
    for(int i=2; i<n+1; ++i) {
      q_fact *= quantum_integer(i);//*= [i]_q
    }
    return q_fact;
  }

/** \brief Return the q-exponential map of a matrix with Ratioanl_mp coefficients, given by the sequence of its non-trivial powers.
 * 
 * The input is the sequence of non-zero powers of x, i.e., pow_x[i] = x^i. 
 * 
 * \f[\operatorname{exp}_q(M) = \displaystyle\sum_{n=0}^{+\infty} \frac{q^{n(n-1)/4}}{[n]_q!} M^n = \displaystyle\sum_{n=0}^{+\infty} \frac{X^{n(n-1)}}{[n]_q!} M^n \f].
 * */
  template<typename Matrix>
  Matrix q_exponential_map(std::vector<Matrix>& pow_x) {
    if(pow_x.empty()) { 
      std::cerr << "Cannot compute the quantum exponential on an empty set of powers.\n"; 
    }
    auto it = pow_x.begin();
    Matrix expq_x = *it;//1/id of appropriate size
    ++it;
    int n=1;
    while(it != pow_x.end()) {
      expq_x += (quantum_monomial(n*(n-1)) / quantum_factorial(n)) * (*it);
      ++it; ++n;
    } 
    return expq_x;
  }

/** \brief Return the q-exponential map of a matrix x.
 * 
 * The input is the sequence of non-zero powers of x, i.e., pow_x[i] = x^i. 
 * 
 * \f[\operatorname{exp}_{q^{-1}}(M) = \displaystyle\sum_{n=0}^{+\infty} \frac{q^{-n(n-1)/4}}{[n]_q!} M^n = \displaystyle\sum_{n=0}^{+\infty} \frac{1}{X^{n(n-1)/2} \times [n]_q!} M^n \f].
 * */
  template<typename Matrix>
  Matrix qinv_exponential_map(std::vector<Matrix>& pow_x) {
    if(pow_x.empty()) { 
      std::cerr << "Cannot compute the quantum exponential on an empty set of powers.\n"; 
    }
    auto it = pow_x.begin();
    Matrix expq_x = *it;//1/id of appropriate size
    ++it;
    int n=1;
    while(it != pow_x.end()) {
      expq_x += (quantum_monomial(n*(1-n)) / quantum_factorial(n)) * (*it);
      ++it; ++n;
    } 
    return expq_x;
  }

// /** \name Methods for the implementation of framed oriented tangle/link invariants derived from an irreducible representation \f$\rho: A \to \operatorname{End}(V)\f$ of a ribbon Hopf algebras \f$A\f$ (with data) into a vector space \f$V\f$, as described in "Ohtsuki - Quantum invariants: A study of knots, 3-manifolds, and their sets" Chapter 4, using the same notations.
//  * @{ */
//   Matrix rho_u(int N, Quantum_group_type Uqg) {
//     switch(Uqg) {
//     case Uqsl2C:
//       break;
//     default: std::cerr << "Invalid quantum group.\n";
//     }
//   }
//   Matrix rho_v(int N, Quantum_group_type Uqg) {
//     switch(Uqg) {
//     case Uqsl2C:
//       break;
//     default: std::cerr << "Invalid quantum group.\n";
//     }
//   }
//   Matrix rho_n(int N, Quantum_group_type Uqg) {
//     switch(Uqg) {
//     case Uqsl2C:
//       break;
//     default: std::cerr << "Invalid quantum group.\n";
//     }
//   }
//   Matrix rho_np(int N, Quantum_group_type Uqg) {
//     switch(Uqg) {
//     case Uqsl2C:
//       break;
//     default: std::cerr << "Invalid quantum group.\n";
//     }
//   }

// /** \brief Return the R-matrix as defined in Ohtsuki (4.30) p.77.*/
//   Matrix h_matrix(int N, Quantum_group_type Uqg) {
//     switch(Uqg) {
//     case Uqsl2C:
//       break;
//     default: std::cerr << "Invalid quantum group.\n";
//     }
//   }
// /** \brief Return the R-matrix as defined in Ohtsuki (4.30) p.77.*/
//   Matrix R_matrix(int N) {

//   }
// /** \brief Return the inverse of the R-matrix for the N-th colored Jones polynomial of a braid.*/
//   Matrix R_inv_matrix(int N) {

//   }


// /* @} */  // end AbelianGroup methods


//   Laurent_poly jones_polynomial(Braid b, int N = 2) {

//   }

// private:
//   Laurent_poly_struct L_;
//   std::map<int, >
};

} //namespace kumquat

#endif // KUMQUAT_JONES_POLYNOMIAL_H_
