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

/**
 * @brief      Class computing and storing data related to quantum group, with a focus to the modular category structure of the representation of the quantum group.
 * 
 * The approach is lazy and keep in memory the useful matrices that have already 
 * been computed.
 * 
 * \implements RibbonCategory.
 * 
 * We follow the conventions of \cite Ohtsuki02 chapter 4.
 */
class Quantum_group_Uqsl2_gen_q {
  /**
   * A multiprecision integer type.
   */
  typedef boost::multiprecision::mpz_int Integer;
  /**
   * Types for rational functions with integer coefficients.
   */
  typedef Rational_function_integral_mp Rational_f;
  /**
   * Morphisms are represented by matrices with rational functions coefficients. For a vector space \f$V\f$ of finite dimension \f$N\f$, we consider a standard basis \f$e_1, \ldots, e_N\f$ in which the matrices are expressed. When considering morphisms \f$V \otimes V \to V \otimes V\f$ we use the basis \f$e_1 \otimes e_1, e_1 \otimes e_2, \ldots, e_1 \otimes e_N, e_2 \otimes e_1, \ldots, e_N \otmies e_N\f$ in that particular order. 
   */
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
/**
 * @brief      Return the n by n identity matrix with rational function coefficients.
 *
 * @param[in]  n     Size of the identity matrix. Must be \f$>0\f$.
 *
 * @return     Return the \f$n\times n\f$ identity matrix.
 */
  Matrix id(int n) {
    Matrix id_n(n,n,Rational_function());
    for(int i=0; i<n; ++i) {
      id_n(i,i) = Rational_f(1);
    }
    return id_n;
  }
  /**
   * @brief      Compute the matrix \f$E\f$, as in \cite Ohtsuki02 chapter 4, p.89.
   *
   * @param[in]  N     The dimension of the irreducible representation of \f$\operatorname{sl}_2(\C)\f$. Must be \f$N > 1\f$.
   *
   * @return     The \f$N\times N\f$ matrix \f$E\f$.
   */
  Matrix E(int N) {
    Matrix e(N,N,Rational_function());
    for(int i=1; i<N; ++i) {
      e(i-1,i) = qd_.quantum_integer(N-i);
    }
    return e;
  }
  /**
   * @brief      Compute the matrix \f$F\f$, as in \cite Ohtsuki02 chapter 4, p.89.
   *
   * @param[in]  N     The dimension of the irreducible representation of \f$\operatorname{sl}_2(\C)\f$. Must be \f$N > 1\f$.
   *
   * @return     The \f$N\times N\f$ matrix \f$F\f$.
   */  
  Matrix F(int N) {
    Matrix f(N,N,Rational_function());
    for(int i=1; i<N; ++i) {
      f(i,i-1) = qd_.quantum_integer(i);
    }
    return f;    
  }
  /**
   * @brief      Compute the matrix \f$K\f$, as in \cite Ohtsuki02 chapter 4, p.89.
   *
   * @param[in]  N     The dimension of the irreducible representation of \f$\operatorname{sl}_2(\C)\f$. Must be \f$N > 1\f$.
   *
   * @return     The \f$N\times N\f$ matrix \f$K\f$.
   */  
  Matrix K(int N) {
    Matrix f(N,N,Rational_function());
    for(int i=0; i<N; ++i) {
      //diag(q^{-(n-1)/2}, q^{-(n-3)/2}, q^{-(n-5)/2}, ...)
      f(i,i) = qd_.quantum_monomial(2* (N - (2*i+1)) );
    }
    return f;    
  }
  /**
   * @brief      Compute the matrix \f$K^{-1}\f$, as in \cite Ohtsuki02 chapter 4, p.89.
   *
   * @param[in]  N     The dimension of the irreducible representation of \f$\operatorname{sl}_2(\C)\f$. Must be \f$N > 1\f$.
   *
   * @return     The \f$N\times N\f$ matrix \f$K^{-1}\f$.
   */  
  Matrix Kinv(int N) {
    Matrix f(N,N,Rational_function());
    for(int i=0; i<N; ++i) {
      //diag(q^{-(n-1)/2}, q^{-(n-3)/2}, q^{-(n-5)/2}, ...)
      f(i,i) = qd_.quantum_monomial((-1)* 2* (N - (2*i+1)) );
    }
    return f;    
  }
/** \brief The matrix P, which is the block N times N matrix:
 * | 0   id|
 * |id   0 |
 * */
  /**
   * @brief      Compute the permutation matrix \f$P\f$, implementing \f$e_i \otimes e_j \to e_j \otimes e_i\f$.
   *
   * @param[in]  N     The dimension of the irreducible representation of \f$\operatorname{sl}_2(\C)\f$. Must be \f$N > 1\f$.
   *
   * @return     The \f$N\times N\f$ matrix \f$P\f$.
   */
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
  /**
   * @brief     Returns the braiding \f$c_{V,V}: V\otimes V \to V \otimes V\f$. 
   *
   * Following \cite Ohtsuki02 chap.4 p.87, this is the \f$R\f$-matrix:
   * \f[ R := P \circ \mathcal{R}, \text{s.t.} 
   * \mathcal{R} := q^{H\otimes H / 4} \times \exp_{q}( (q^{1/2} - q^{-1/2}) E \otimes F ).\f]
   * 
   * @param[in]  N     The dimension of the irreducible representation of \f$\operatorname{sl}_2(\C)\f$. Must be \f$N > 1\f$.
   *
   * @return     The \f$N\times N\f$ matrix \f$c_{V,V}\f$.
   */
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
  /**
   * @brief     Returns the inverse braiding \f$c^{-1}_{V,V}: V\otimes V \to V \otimes V\f$. 
   *
   * Following \cite Ohtsuki02 chap.4, this is the inverse of the \f$R\f$-matrix:
   * \f[ R^{-1} := \mathcal{R^{-1}} \circ P, \text{s.t.} 
   * \mathcal{R^{-1}} := \exp_{q^{-1}}( (q^{-1/2} - q^{1/2}) E \otimes F )\times q^{-H\otimes H / 4}\f]
   *
   * @param[in]  N     The dimension of the irreducible representation of \f$\operatorname{sl}_2(\C)\f$. Must be \f$N > 1\f$.
   *
   * @return     The \f$N\times N\f$ matrix \f$c^{-1}_{V,V}\f$.
   */
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







private:
  //p.88
  Morphism u(int N) {

    INVESTIGATE:
    /* (FK^-1)^n  =  q^{n(n-1)/2}F^n K^-n
     * the power does not act pointwise as usual. p.337
     *
     *
     *
     *
     *
     *
     */

    //exp_q part
    //prepare powers all_pow[k] =  q^{2n(n-1)/k} (q^{-1/2}-q^{1/2})^k (FK^{-1}E)^k
    std::vector< Morphism > all_pow; all_pow.reserve(N-1);
    all_pow.push_back(id(N*N));//N^2 x N^2 id
    auto E = E(); auto F = F(); auto Kinv = Kinv();

    F *= (-1) * qd_.quantum_half();// (q^{-1/2}-q^{1/2})F
    F *= Kinv;
    F *= E; //(q^{-1/2}-q^{1/2})(FK^{-1}E)
    all_pow.push_back(F);
    auto pow_FKE = F;//==(q^{-1/2}-q^{1/2})(FK^{-1}E)
    for(int i=2; i<N; ++i) {
      pow_FKE *= F; //[(q^{-1/2}-q^{1/2})(FK^{-1}E)]^{i}
      all_pow.push_back(qd_.quantum_monomial(2*i*(i-1)) * pow_FKE);
    }
    auto u_mat = qd_.q_exponential_map(all_pow);

    //multiply on the left by q^{-H^2/4} = diag(q^{-(N-1)^2}, q^{-(N-3)^2}, ...)
    // row_i(u_mat) *= q^{-(N - (2i+1))^2} / 4 = X^{-( N - (2i+1) )^2}
    for(int i=0; i<N; ++i) {//row i
      for(int j=0; j<N; ++j) {//col j
        u_mat(i,j) *= qd_.quantum_monomial( (-1)* ( N- 2*i -1 ) * ( N- 2*i -1 ));
      }
    }
    return u_mat;
  }
  //p.88
  Morphism v(int N) {
    return (Kinv()*u());
  }

  //p.77
  Morphism h(int N) {
    return K();
  }




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