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

#ifndef KUMQUAT_QUANTUM_GROUP_H_ 
#define KUMQUAT_QUANTUM_GROUP_H_

#include <kumquat/Vector_space.h>
#include <kumquat/Dense_matrix.h>
#include <kumquat/Rational_function_integral_mp.h>

namespace kumquat {

// enum Quantum_group_type {Uq_sl2C_generic_q, Uq_sl2C_rootofunit_q};

/** \brief Class computing and storing data related to quantum group, with a focus to the modular category structure of the representation of the quantum group.
 * 
 * The approach is lazy and keep in memory the useful matrices that have already 
 * been computed.
 * 
 * Is model of concept RibbonCategory
 **/
// template<int N>
class Quantum_group_Uqsl2_gen_q {

  typedef boost::multiprecision::mpz_int Integer;
  typedef Rational_function_integral_mp Rational_f;
/** \brief Morphisms are represented by matrices with rational functions coefficients.*/
  typedef Dense_matrix<Rational_f> Morphism;

/** \brief A handle type to designate an object in the category. Objects are finite dimensional vector spaces of a field (generally \$f\mathbb{C}\f$).*/
  // typedef unspecified Object_handle;
/** \brief A handle type to designate a morphism in the category. Morphism are matrices of dimension compatible with the vectors.*/
  // typedef unspecified Morphism_handle;

  Quantum_group_Uqsl2_gen_q() {
    qd_();
    ;
  }

/** \brief Return the n by n identity matrix with rational function coefficients.*/
  Matrix id(int n) {
    Matrix id_n(n,n,Rational_function());
    for(int i=0; i<n; ++i) {
      id_n(i,i) = Rational_f(1);
    }
    return id_n;
  }
/** \brief The matrix E.*/
  Matrix E(int N) {
    Matrix e(N,N,Rational_function());
    for(int i=1; i<N; ++i) {
      e(i-1,i) = qd_.quantum_integer(N-i);
    }
    return e;
  }
/** \brief The matrix E.*/
  Matrix F(int N) {
    Matrix f(N,N,Rational_function());
    for(int i=1; i<N; ++i) {
      f(i,i-1) = qd_.quantum_integer(i);
    }
    return f;    
  }
/** \brief The matrix P, which is the block N times N matrix:
 * | 0   id|
 * |id   0 |
 * */
  Matrix P(int N) {
    Matrix p(N*N,N*N,Rational_function());
    for(int i=0; i<N; ++i) {
      for(int j=0; j<N; ++j) {
        //e_i \otimes e_j -> e_j \otimes e_i
        P(i*N+j,j*N+i) = Rational_f(1);
      }
    }
    return P;
  }
/** \brief Returns the braiding \f$c_{V,W}: V\otimes W \to W \otimes V\f$.
 * 
 * Input the handle for objects \f$V\f$ and \f$W\f$.*/  
  // Morphism_handle braiding(Object_handle v_h, Object_handle w_h)  
  Matrix braiding(int N) {
    //only one simple object <-> irreducible rep of dim N

    //compute exp_q( (q^{0.5}-q^{-0.5}) E \otimes F )
    //powers of (E \otimes F)^k = E^k \otimes F^k;
    std::vector< Morphism > pow_EoF; all_pow_EoF.reserve(N-1);
    all_pow_EoF.push_back(id(N*N));//N^2 x N^2 id
    auto E = E(); auto F = F();
    E *= qd_.quantum_half();// (q^{1/2}-q^{-1/2}E)
    all_pow_EoF.push_back(E.rtensor(F));
    auto pow_E = E;    auto pow_F = F;
    for(int i=2; i<N; ++i) {
      pow_E *= E; pow_F *= F;//E^i, F^i
      all_pow_EoF.push_back(pow_E.rtensor(pow_F));
    }
    auto R = qd_.q_exponential_map(all_pow_EoF);//==exp_q( (q^{0.5}-q^{-0.5}) E \otimes F )

    //compute the part q^{1/4 * H \otimes H} = exp(h/4 * H \otimes H) = X^{H \otimes H}, with H = diag(n-1, n-3, n-5, ..., -(n-1))
    std::vector<Rational_f> diag_X_to_HoH(N*N);
    for(int i=0; i<N; ++i) {
      for(int j=0; j<N; ++j) {
        //diagonal coefficients of X^{H \otimes H}
        diag_X_to_HoH(N*i+j) = Q_data_.quantum_monomial((N-(2*i)-1)*(N-(2*j)-1));
      }
    }
    //multiply R on the left by X^{H \otimes H} 
    for(int i=0; i<N*N; ++i) {
      for(int j=0; j<N*N; ++j) {
        R(i,j) *= diag_X_to_HoH[i]; 
      }
    }
    R.ltimes_equal(P());
    return R;
  }
/** \brief Returns the inverse braiding \f$c^{-1}_{V,W}: W \otimes V \to V \otimes W\f$.
 * 
 * 
 * Return the R matrix such that:
 * \f[
 * R^{-1} = \mathcal{R^{-1}} \circ P, \mathcal{R^{-1}} = \exp_{q^{-1}}( (q^{-1/2} - q^{1/2}) E \otimes F )\times q^{-H\otimes H / 4}
 * \f]
 * */
  // Morphism_handle braiding_inv(Object_handle v_h, Object_handle w_h) {
  Morphism braiding_inv(int N) {
    //only one simple object <-> irreducible rep of dim N

    //compute exp_q( (q^{0.5}-q^{-0.5}) E \otimes F )
    //powers of (E \otimes F)^k = E^k \otimes F^k;
    std::vector< Morphism > pow_EoF; all_pow_EoF.reserve(N-1);
    all_pow_EoF.push_back(id(N*N));//N^2 x N^2 id
    auto E = E(); auto F = F();
    E *= qd_.quantum_half();// (q^{1/2}-q^{-1/2})E
    E *= -1;// (q^{-1/2}-q^{1/2})E
    all_pow_EoF.push_back(E.rtensor(F));//(q^{-1/2}-q^{1/2})E \otimes F)
    auto pow_E = E;    auto pow_F = F;
    for(int i=2; i<N; ++i) {
      pow_E *= E; pow_F *= F;//(q^{-1/2}-q^{1/2})^iE^i, F^i
      all_pow_EoF.push_back(pow_E.rtensor(pow_F));
    }
    auto Rinv = qd_.qinv_exponential_map(all_pow_EoF);//==exp_{q^{-1}}( (q^{-1/2}-q^{1/2}) E \otimes F )

    //compute the part q^{- 1/4 * H \otimes H} = exp(-h/4 * H \otimes H) = X^{-H \otimes H}, with H = diag(n-1, n-3, n-5, ..., -(n-1))
    std::vector<Rational_f> diag_X_to_HoH(N*N);
    for(int i=0; i<N; ++i) {
      for(int j=0; j<N; ++j) {
        //diagonal coefficients of X^{-H \otimes H}
        diag_X_to_HoH(N*i+j) = Q_data_.quantum_monomial((-1)*(N-(2*i)-1)*(N-(2*j)-1));
      }
    }
    //multiply R^{-1} on the right by X^{H \otimes H} 
    for(int i=0; i<N*N; ++i) {
      for(int j=0; j<N*N; ++j) {
        Rinv(i,j) *= diag_X_to_HoH[i]; 
      }
    }
    Rinv.rtimes_equal(P());
    return Rinv;
  }

/** \brief Returns the twist morphism \f$\theta_V: V \to V\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism twist(int N) {
  }

/** \brief Returns the inverse of the twist morphism \f$\theta^{-1}_V: V \to V\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism twist_inv(int N) {
    
  }
/** \brief Return the dual of an object \fV^*\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Object_handle dual(Object_handle v_h);
/** \brief Return the pairing morphism \f$d_V: V^* \otimes V \to \mathbbm{1}\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism_handle pairing(Object_handle v_h);
/** \brief Return the copairing morphism \f$b_V: \mathbbm{1} \to V \otimes V^* \to \f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism_handle copairing(Object_handle v_h);
/** \brief Return the dimension of an object.
 *
 * Input the handle for object \f$V\f$, return an element of \f$\operatorname{End}(\mathbbm{1})\f$.*/
  Base_element dim(Object_handle v_h);
/** Return the trace of a morphism.*/
  Base_element trace(const Morphism& phi);


// /** \brief Returns the twist morphism \f$\theta_V: V \to V\f$.
//  * 
//  * Input the handle for object \f$V\f$.*/
//   Morphism_handle twist(Object_handle v_h);
// /** \brief Returns the inverse of the twist morphism \f$\theta^{-1}_V: V \to V\f$.
//  * 
//  * Input the handle for object \f$V\f$.*/
//   Morphism_handle twist_inv(Object_handle v_h);
// /** \brief Return the dual of an object \fV^*\f$.
//  * 
//  * Input the handle for object \f$V\f$.*/
//   Object_handle dual(Object_handle v_h);
// /** \brief Return the pairing morphism \f$d_V: V^* \otimes V \to \mathbbm{1}\f$.
//  * 
//  * Input the handle for object \f$V\f$.*/
//   Morphism_handle pairing(Object_handle v_h);
// /** \brief Return the copairing morphism \f$b_V: \mathbbm{1} \to V \otimes V^* \to \f$.
//  * 
//  * Input the handle for object \f$V\f$.*/
//   Morphism_handle copairing(Object_handle v_h);
// /** \brief Return the dimension of an object.
//  *
//  * Input the handle for object \f$V\f$, return an element of \f$\operatorname{End}(\mathbbm{1})\f$.*/
//   Base_element dim(Object_handle v_h);
// /** Return the trace of a morphism.*/
//   Base_element trace(const Morphism& phi);


  void compute_H() {

  }



private:
  Morphism H_;
  Quantum_data qd_;

  std::map<int,Matrix> pow_R_;
  std::map<int,Matrix> pow_Rinv_;
  std::map<int,Matrix> pow_twist_;
  std::map<int,Matrix> pow_twistinv_;
  std::map<int,Matrix> pow_pairing_;
  std::map<int,Matrix> pow_copairing_;//and duals?

};

} // namespace kumquat

#endif // END KUMQUAT_QUANTUM_GROUP_H_ 