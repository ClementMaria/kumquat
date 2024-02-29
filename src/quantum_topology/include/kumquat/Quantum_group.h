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

#include <kumquat/Quantum_data.h>
#include <kumquat/Dense_matrix.h>
#include <kumquat/Rational_function_integral_mp.h>
#include <kumquat/Braid.h>
#include <kumquat/Plat_braid.h>
#include <kumquat/Rational_function.h>

namespace kumquat {

/**
 * @brief      Class computing and storing data related to quantum group, with a focus to the modular category structure of the representation of the quantum group.
 * 
 * The approach is lazy and keep in memory the useful matrices that have already 
 * been computed.
 * 
 * \implements RibbonCategory.
 *
 * @tparam     N     The dimension of the irreducible representation of \f$U_q(\mathoperator{sl}_2(\mathbb{C}))\f$.
 * 
 * We follow the conventions of \cite Ohtsuki02 chapter 4.
 */
template<int N>
class Quantum_group_Uqsl2_gen_q {
public:
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
  typedef Dense_matrix<Rational_function> Morphism;

/** \brief A handle type to designate an object in the category. Objects are finite dimensional vector spaces of a field (generally \$f\mathbb{C}\f$).*/
  // typedef unspecified Object_handle;
/** \brief A handle type to designate a morphism in the category. Morphism are matrices of dimension compatible with the vectors.*/
  // typedef unspecified Morphism_handle;

  Quantum_group_Uqsl2_gen_q(int num_strands, int max_twist) {
    // h_tensor_(0,0,Rational_function());
    h_tensor_ = h_morphism(num_strands);

    pow_R_.reserve(max_twist+1); 
    pow_R_.push_back(id_morphism(N*N)); pow_R_.push_back(braiding());
    pow_Rinv_.reserve(max_twist+1);
    pow_Rinv_.push_back(id_morphism(N*N)); pow_Rinv_.push_back(braiding_inv());
    for(int i=2; i<= max_twist; ++i) {
      
      // std::cout << " ____________________ Multiplication R * R^" << i-1 << "\n";
      // std::cout << "R == \n";
      // std::cout << pow_R_[1] << "\nand R^" << i-1 << "==\n";
      // std::cout << pow_R_[i-1] << "\n\n";
      // std::cout << "RESULTAT R^" << i << " == \n";
      pow_R_.push_back(pow_R_[1]*pow_R_[i-1]);
      // std::cout << pow_R_[i] << "\n";
      // std::cout << "________________________________________________________";
      pow_Rinv_.push_back(pow_Rinv_[1]*pow_Rinv_[i-1]);
    }
  }

  void display() {
    std::cout << "***** R matrices:--------\n";
    for(size_t i=0; i<pow_R_.size(); ++i) {
      std::cout << "   R^" << i << " = \n";
      std::cout << pow_R_[i] << "\n\n";
    }
    for(size_t i=0; i<pow_Rinv_.size(); ++i) {
      std::cout << "   R^(-" << i << ") = \n";
      std::cout << pow_Rinv_[i] << "\n\n";
    }

    for(size_t i=0; i<pow_R_.size(); ++i) {
      auto m1 = pow_R_[i];
      auto m2 = pow_Rinv_[i];
      m1 *= m2;
      std::cout << "   R^" << i << "*R^(-" << i << ") = \n";
      std::cout << m1 << "\n\n";
    }



  
    std::cout << "***** h tensor:-------\n";
    std::cout << h_tensor_ << "\n\n";
  }

private:
  /**
   * @brief      Return the n by n identity matrix with rational function coefficients.
   *
   * @param[in]  n     Size of the identity matrix. Must be \f$>0\f$.
   *
   * @return     Return the \f$n\times n\f$ identity matrix.
   */
  Morphism id_morphism(int n) {
    Morphism id_n(n,n,Rational_function());
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
  Morphism E_morphism() {
    Morphism e(N,N,Rational_function());
    for(int i=1; i<N; ++i) {
      e(i-1,i) = q_data_.quantum_integer(N-i);
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
  Morphism F_morphism() {
    Morphism f(N,N,Rational_function());
    for(int i=1; i<N; ++i) {
      f(i,i-1) = q_data_.quantum_integer(i);
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
  Morphism K_morphism() {
    Morphism f(N,N,Rational_function());
    for(int i=0; i<N; ++i) {
      //diag(q^{-(n-1)/2}, q^{-(n-3)/2}, q^{-(n-5)/2}, ...)
      f(i,i) = q_data_.quantum_monomial(2* (N - (2*i+1)) );
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
  Morphism Kinv_morphism() {
    Morphism f(N,N,Rational_function());
    for(int i=0; i<N; ++i) {
      //diag(q^{-(n-1)/2}, q^{-(n-3)/2}, q^{-(n-5)/2}, ...)
      f(i,i) = q_data_.quantum_monomial((-1)* 2* (N - (2*i+1)) );
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
  Morphism P_morphism() {
    Morphism P(N*N,N*N,Rational_function());
    for(int i=0; i<N; ++i) {
      for(int j=0; j<N; ++j) {
        //e_i \otimes e_j -> e_j \otimes e_i
        P(i*N+j,j*N+i) = Rational_f(1);
      }
    }
    return P;
  }

public:
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
  Morphism braiding() {
    //only one simple object <-> irreducible rep of dim N
    //compute exp_q( (q^{0.5}-q^{-0.5}) E \otimes F )
    //powers of (E \otimes F)^k = E^k \otimes F^k;
    std::vector< Morphism > all_pow_EoF; all_pow_EoF.reserve(N-1);
    all_pow_EoF.push_back(id_morphism(N*N));//N^2 x N^2 id
    auto E = E_morphism(); auto F = F_morphism();
    E *= q_data_.quantum_half();// (q^{1/2}-q^{-1/2}E)
    all_pow_EoF.push_back(E.rtensor(F));
    auto pow_E = E;    auto pow_F = F;
    for(int i=2; i<N; ++i) {
      pow_E *= E; pow_F *= F;//E^i, F^i
      all_pow_EoF.push_back(pow_E.rtensor(pow_F));
    }
    auto R = q_data_.q_exponential_map(all_pow_EoF);//==exp_q( (q^{0.5}-q^{-0.5}) E \otimes F )

    //compute the part q^{1/4 * H \otimes H} = exp(h/4 * H \otimes H) = X^{H \otimes H}, with H = diag(n-1, n-3, n-5, ..., -(n-1))
    std::vector<Rational_f> diag_X_to_HoH(N*N);
    for(int i=0; i<N; ++i) {
      for(int j=0; j<N; ++j) {
        //diagonal coefficients of X^{H \otimes H}
        diag_X_to_HoH[N*i+j] = q_data_.quantum_monomial((N-(2*i)-1)*(N-(2*j)-1));
      }
    }
    //multiply R on the left by X^{H \otimes H} 
    for(int i=0; i<N*N; ++i) {
      for(int j=0; j<N*N; ++j) {
        R(i,j) *= diag_X_to_HoH[i]; 
      }
    }
    auto P = P_morphism();
    R.ltimes_equal(P);
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
  Morphism braiding_inv() {
    // auto it = pow_braiding_.find(-1);
    // if(it != pow_braiding_.end()) { return *it; }
    //only one simple object <-> irreducible rep of dim N

    //compute exp_q( (q^{0.5}-q^{-0.5}) E \otimes F )
    //powers of (E \otimes F)^k = E^k \otimes F^k;
    std::vector< Morphism > all_pow_EoF; all_pow_EoF.reserve(N-1);
    all_pow_EoF.push_back(id_morphism(N*N));//N^2 x N^2 id
    auto E = E_morphism(); auto F = F_morphism();
    E *= (-1)*q_data_.quantum_half();// (q^{1/2}-q^{-1/2})E
    // E *= (int)(-1);// (q^{-1/2}-q^{1/2})E
    all_pow_EoF.push_back(E.rtensor(F));//(q^{-1/2}-q^{1/2})E \otimes F)
    auto pow_E = E;    auto pow_F = F;
    for(int i=2; i<N; ++i) {
      pow_E *= E; pow_F *= F;//(q^{-1/2}-q^{1/2})^iE^i, F^i
      all_pow_EoF.push_back(pow_E.rtensor(pow_F));
    }
    auto Rinv = q_data_.qinv_exponential_map(all_pow_EoF);//==exp_{q^{-1}}( (q^{-1/2}-q^{1/2}) E \otimes F )

    //compute the part q^{- 1/4 * H \otimes H} = exp(-h/4 * H \otimes H) = X^{-H \otimes H}, with H = diag(n-1, n-3, n-5, ..., -(n-1))
    std::vector<Rational_f> diag_X_to_HoH(N*N);
    for(int i=0; i<N; ++i) {
      for(int j=0; j<N; ++j) {
        //diagonal coefficients of X^{-H \otimes H}
        diag_X_to_HoH[N*i+j] = q_data_.quantum_monomial((-1)*(N-(2*i)-1)*(N-(2*j)-1));
      }
    }
    //multiply R^{-1} on the right by X^{H \otimes H} 
    for(int i=0; i<N*N; ++i) {
      for(int j=0; j<N*N; ++j) {
        Rinv(i,j) *= diag_X_to_HoH[i]; 
      }
    }
    auto P = P_morphism();
    Rinv.rtimes_equal(P);
    // pow_braiding_inv_[1] = Rinv;//remember the computation
    return Rinv;
  }

/** \brief Returns the twist morphism \f$\theta_V: V \to V\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  // Morphism twist() {}
/** \brief Returns the inverse of the twist morphism \f$\theta^{-1}_V: V \to V\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  // Morphism twist_inv() {}
/** \brief Return the dual of an object \fV^*\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  // Object_handle dual(Object_handle v_h);
/** \brief Return the pairing morphism \f$d_V: V^* \otimes V \to \mathbbm{1}\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  // Morphism_handle pairing(Object_handle v_h);
/** \brief Return the copairing morphism \f$b_V: \mathbbm{1} \to V \otimes V^* \to \f$.
 * 
 * Input the handle for object \f$V\f$.*/
  // Morphism_handle copairing(Object_handle v_h);
/** \brief Return the dimension of an object.
 *
 * Input the handle for object \f$V\f$, return an element of \f$\operatorname{End}(\mathbbm{1})\f$.*/
  // Rational_f dim(Object_handle v_h);
/** Return the trace of a morphism.*/
  Rational_f trace(Morphism& phi) {
    return phi.trace();
  }

  //p.77
  //return h^{\otimes n}
  Morphism h_morphism(int tensor_n) {
    // auto it = tensors_h_.find(tensor_n);
    // if(it != tensors_h_.end()) { return *it; }
    Morphism h = K_morphism();
    Morphism tens_h = h;
    for(int i=2; i<= tensor_n; ++i) {
      tens_h.rtensor_equal(h);
    }
    // tensors_h_[tensor_n] = tens_h;
    return tens_h;
  }

/**
 * @brief      Compute the quantum invariant on an input braid, associated with the ribbon category.
 *
 * @param[in]  b     A compact description of a braid.
 *
 * @return     The quantum invariant as a rational function.
 */
  Rational_f quantum_invariant(Braid& b) {

    if(b.braid().empty()) { return 0; }

    int num_strands = b.num_strands();
    auto it = b.braid().begin();
    Morphism tau = id_R_id(it->first,it->second,num_strands);
    ++it;

    std::cout << "+++++++++++++++++++++++++++++++++\n\n";
    std::cout << tau << "\n";
    std::cout << "+++++++++++++++++++++++++++++++++\n\n";

    for(; it!=b.braid().end(); ++it) {
      std::cout << "(" << it->first << "," << it->second << "," << num_strands << ")\n";
      tau.ltimes_equal(id_R_id(it->first,it->second,num_strands));
    }
    tau.ltimes_equal(h_tensor_);
    return trace(tau);
  }


/**
 * @brief      Compute the quantum invariant on an input plat braid, and take the closing..., associated with the ribbon category.
 *
 * @param[in]  b     A description of a plat braid.
 *
 * @return     The quantum invariant as a rational function.
 */
  Rational_f quantum_invariant(Plat_braid& b) {

    if(b.braid().empty()) { return 0; }

    int num_strands = b.num_strands();
    auto it = b.braid().begin();
    Morphism tau = id_R_id(it->begin(),it->end(),num_strands);
    ++it;

    // std::cout << "Size tau = " << tau.num_rows() << " x " << tau.num_columns() << "\n";


    for(; it!=b.braid().end(); ++it) {

      // for(auto pp : *it) {
      //   std::cout << "[" << pp.first << "," << pp.second <<"] ";
      // }
      // std::cout << std::endl;

      tau.ltimes_equal(id_R_id(it->begin(),it->end(),num_strands));
    
      // std::cout << "Size tau = " << tau.num_rows() << " x " << tau.num_columns() << "\n";

    }

    // std::cout << "B\n";

    // std::cout << "Size tau = " << tau.num_rows() << " x " << tau.num_columns() << "     size h = " << h_tensor_.num_rows() << " x " << h_tensor_.num_columns() <<"\n";

    tau.ltimes_equal(h_tensor_);
    
    // std::cout << "C\n";

    return trace(tau);
  }
  /**
   * @brief      Compute and return the matrix for the morphism corresponding to a collection of independent parallel twist regions.
   *
   * If the iterators point towards the pairs \f$(i_1,n_1), \ldots, (i_k,n_k)\f$, with \f$n_j > 0\f$ and \f$i_j + 2 \leq i_{j+1}\f$ for all appropriate \f$j\f$, the output morphism is:
   * \f[
   *      \operatorname{id}_N^{\otimes (|i_1|-1)} \otimes R^{\operatorname{sign}(i_1) n_1} \otimes \operatorname{id}_N^{\otimes (|i_2|-|i_1|-2)} \otimes R^{\operatorname{sign}(i_2) n_2} \otimes \ldots
   * \f]
   * where \f$\operatorname{id}_N\f$ is the N by N identity. 
   *
   * @param[in]  beg            The beg
   * @param[in]  end            The end
   * @param[in]  num_strands    The number of strands in the braid.
   *
   * @tparam     IteratorPairs  Iterator type on a set of ordered pairs (i,n), n>0, representing twist regions as in the Plat_braid structure.
   *
   * @return     The morphism.
   */
  template<typename IteratorPairs >
  Morphism id_R_id(IteratorPairs beg, IteratorPairs end, int num_strands) {
    if(beg == end) { return id_morphism(std::pow(N,num_strands)); }
    Morphism M;//(0,0,Rational_function());
    auto curr_idx = beg->first;
    //initialize to \operatorname{id}_N^{\otimes (i_1-1)} \otimes R^{n_1}
    if(beg->second < 0) {//n_1 < 0 => negative twists
      M = (pow_Rinv_[std::abs(beg->second)]).ltensor(id_morphism(std::pow(N,curr_idx)));
    }
    else {//n_1 > 0 => positive twists
      M = (pow_R_[beg->second]).ltensor(id_morphism(std::pow(N,curr_idx)));
    }
    ++beg;
    auto prev_idx = curr_idx;
    while(beg != end) {
      curr_idx = beg->first;
      if(beg->second < 0) {//n_i < 0
        M.rtensor_equal( id_morphism(std::pow(N,curr_idx-prev_idx-2)).rtensor(pow_Rinv_[std::abs(beg->second)]) );
      }
      else {//n_i > 0
        M.rtensor_equal( id_morphism(std::pow(N,curr_idx-prev_idx-2)).rtensor(pow_R_[beg->second]) );
      }
      prev_idx = curr_idx;
      ++beg;
    }
    //add the last \otimes id_N^{num_strands-|i_k|-1}
    M.rtensor_equal(id_morphism(std::pow(N,num_strands-curr_idx-2)));
    return M;
  }

/**
 * @brief      Compute and return the matrix for the morphism 
 * \f[
 *      \operatorname{id}_N^{\otimes (i-1)} \otimes R^{k} \otimes \operatorname{id}_N^{\otimes (n-i-1)},
 * \f]
 * where \f$\operatorname{id}_N\f$ is the N by N identity. 
 *
 * @param[in]  strand_idx    The index \f$i\f$ of the strand, to indicate a crossing between strand |i|>0 and strand |i|+1. If i>0 the crossing is positive, and if i<0 the crossing is negative.
 * @param[in]  twist_length  The number of consecutive crossings in the twist region.
 * @param[in]  num_strands   The total number of strands in the braid.
 *
 * @return     The morphism.
 */
  Morphism id_R_id(int strand_idx, int twist_length, int num_strands) {
    if(twist_length < 0) {
      return (pow_Rinv_[std::abs(twist_length)].ltensor(id_morphism(std::pow(N,strand_idx))) ).rtensor(id_morphism(std::pow(N,num_strands-strand_idx-2)));
    }
    return (pow_R_[twist_length].ltensor(id_morphism(std::pow(N,strand_idx))) ).rtensor(id_morphism(std::pow(N,num_strands-strand_idx-2)));
  }

/**
 * @brief      Return the quantum invariant of the unknot.
 *
 * @return     The quantum invariant of the unknot, as a rational function.
 */
  Rational_f quantum_invariant_unknot() {
    Rational_f sum(0);
    for(int i=0; i<N; ++i) {
      sum += q_data_.quantum_monomial(2*(N-1 -2*i));
    }
    return sum;
  }

private:
  //p.88 and p.337
  //
  Morphism u() {
    //exp_q part
    //prepare powers all_pow[n] =  (q^{-1/2}-q^{1/2})^n (FK^{-1})^n E^n
    std::vector< Morphism > all_pow; all_pow.reserve(N-1);
    all_pow.push_back(id_morphism(N));//N x N id for n=0
    auto E = E_morphism(); auto FKinv = F_morphism()*Kinv_morphism();

    FKinv *= (-1) * q_data_.quantum_half();// (q^{-1/2}-q^{1/2})FK^{-1}
    all_pow.push_back(E*FKinv);//n=1
    auto pow_FKinv = FKinv;//==(q^{-1/2}-q^{1/2})(FK^{-1}E)
    auto pow_E = E;
    for(int i=2; i<N; ++i) {
      pow_FKinv *= FKinv; //(q^{-1/2}-q^{1/2})^i (FK^{-1})^i
      pow_E *= E;//E^i
      all_pow.push_back(pow_FKinv * pow_E);//pow[n] = (q^{-1/2}-q^{1/2})^i (FK^{-1})^i E^i
    }
    //sum_n  q^{n(n-1)/4}/[n]!  (q^{-1/2}-q^{1/2})^n (FK^{-1})^n E^n
    auto u_mat = q_data_.q_exponential_map(all_pow);

    //multiply on the left by q^{-H^2/4} = diag(q^{-(N-1)^2}, q^{-(N-3)^2}, ...)
    //row_i(u_mat) *= q^{-(N - (2i+1))^2} / 4 = X^{-( N - (2i+1) )^2}
    for(int i=0; i<N; ++i) {//row i
      for(int j=0; j<N; ++j) {//col j
        u_mat(i,j) *= q_data_.quantum_monomial( (-1)* ( N- 2*i -1 ) * ( N- 2*i -1 ));
      }
    }
    return u_mat;
  }
  //p.88
  Morphism v() {
    return (Kinv_morphism()*u());
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
//   Rational_f dim(Object_handle v_h);
// /** Return the trace of a morphism.*/
//   Rational_f trace(const Morphism& phi);


private:
  Quantum_data<Rational_function_integral_mp> q_data_;
/**
 * A certain tensor of h, fixed by the number of strands considered.
 */
  Morphism h_tensor_;

  std::vector<Morphism> pow_R_;
  std::vector<Morphism> pow_Rinv_;

/**
 * The h morphism, see \cite Ohtsuki02 Chapter 4 p. 89 (==K)
 */
  // Morphism h_;
  /**
   * tensors_h_[n] = h^{\otimes n}
   */
  // std::map<int,Morphism> tensors_h_;
/**
 * pow_braiding_[n] = R^n. Computed lazily.
 */
  // std::map<int,Morphism> pow_R_;
  // std::map<int,Morphism> pow_Rinv_;
  // std::map<int,Morphism> pow_twist_;
  // std::map<int,Morphism> pow_twistinv_;
  // std::map<int,Morphism> pow_pairing_;
  // std::map<int,Morphism> pow_copairing_;//and duals?

};

} // namespace kumquat

#endif // END KUMQUAT_QUANTUM_GROUP_H_ 