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

#include <kumquat/Laurent_polynomial_mp.h>
#include <kumquat/Dense_matrix.h>
#include <kumquat/Braid.h>
#include <kumquat/Plat_braid.h>
#include <kumquat/Ring_over.h>

namespace kumquat {

/**
 * @brief      
 * 
 * \implements RibbonCategory.
 */
class Jones_polynomial {
public:
  /**
   * Types for rational functions with integer coefficients.
   */
  typedef Laurent_polynomial_mp Laurent_p;
  /**
   * A multiprecision integer type.
   */
  typedef Laurent_p::Coefficient Coeff_p; 
  /**
   * Morphisms are represented by matrices with rational functions coefficients. For a vector space \f$V\f$ of finite dimension \f$N\f$, we consider a standard basis \f$e_1, \ldots, e_N\f$ in which the matrices are expressed. When considering morphisms \f$V \otimes V \to V \otimes V\f$ we use the basis \f$e_1 \otimes e_1, e_1 \otimes e_2, \ldots, e_1 \otimes e_N, e_2 \otimes e_1, \ldots, e_N \otmies e_N\f$ in that particular order. 
   */
  typedef Dense_matrix< Ring_over<Laurent_p> > Morphism;

/** \brief A handle type to designate an object in the category. Objects are finite dimensional vector spaces of a field (generally \$f\mathbb{C}\f$).*/
  // typedef unspecified Object_handle;
/** \brief A handle type to designate a morphism in the category. Morphism are matrices of dimension compatible with the vectors.*/
  // typedef unspecified Morphism_handle;

  Jones_polynomial(int num_strands, int max_twist) {
    N=2;
    // h_tensor_(0,0,Laurent_punction());
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

    std::cout << "--- R = \n" << pow_R_[1] << "\n";
    std::cout << "--- R^-1 = \n" << pow_Rinv_[1] << "\n";


    // display();
  }

  void display() {
    // std::cout << "***** R matrices:--------\n";
    // for(size_t i=0; i<pow_R_.size(); ++i) {
    //   std::cout << "   R^" << i << " = \n";
    //   std::cout << pow_R_[i] << "\n\n";
    // }
    // for(size_t i=0; i<pow_Rinv_.size(); ++i) {
    //   std::cout << "   R^(-" << i << ") = \n";
    //   std::cout << pow_Rinv_[i] << "\n\n";
    // }

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
    Morphism id_n(n,n,Ring_over<Laurent_p>());
    for(int i=0; i<n; ++i) {
      id_n(i,i) = Laurent_p(0,1);
    }
    return id_n;
  }

public:
  /**
   * @brief     Ohtsuki p.30, t^0.5 = X. 
   */
  Morphism braiding() {
    Dense_matrix R(4,4,Ring_over<Laurent_p>());
    R(0,0) = Laurent_p(1,1);
    R(2,1) = Laurent_p(2,1);
    R(1,2) = Laurent_p(2,1);
    R(2,2) = (Laurent_p(1,1) + Laurent_p(3,-1));
    R(3,3) = Laurent_p(1,1);
    return R;
  }
  /**
   * @brief     Ohtsuki p.32
   */
  Morphism braiding_inv() {
    Morphism Rinv(4,4,Ring_over<Laurent_p>());
    Rinv(0,0) = Laurent_p(-1,1);
    Rinv(1,1) = (Laurent_p(-3,-1) + Laurent_p(-1,1));
    Rinv(1,2) = Laurent_p(-2,1);
    Rinv(2,1) = Laurent_p(-2,1);
    Rinv(3,3) = Laurent_p(-1,1);
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
  // Laurent_p dim(Object_handle v_h);
/** Return the trace of a morphism.*/
  Laurent_p trace(Morphism& phi) {
    return phi.trace();
  }

  //p.77
  //return h^{\otimes n}
  Morphism h_morphism(int tensor_n) {
    // auto it = tensors_h_.find(tensor_n);
    // if(it != tensors_h_.end()) { return *it; }
    Morphism h(2,2,Ring_over<Laurent_p>());
    h(0,0) = Laurent_p(-1,1);
    h(1,1) = Laurent_p(1,1);
    Morphism tens_h = h;
    for(int i=2; i<= tensor_n; ++i) {
      tens_h.rtensor_equal(h);
    }
    return tens_h;
  }

/**
 * @brief      Compute the quantum invariant on an input braid, associated with the ribbon category.
 *
 * @param[in]  b     A compact description of a braid.
 *
 * @return     The quantum invariant as a rational function.
 */
  Laurent_p quantum_invariant(Braid& b) {

    if(b.braid().empty()) { return Laurent_p(); }

    int num_strands = b.num_strands();
    auto it = b.braid().begin();
    Morphism tau = id_R_id(it->first,it->second,num_strands);
    ++it;

    // std::cout << "+++++++++++++++++++++++++++++++++\n\n";
    // std::cout << tau << "\n";
    // std::cout << "+++++++++++++++++++++++++++++++++\n\n";

    for(; it!=b.braid().end(); ++it) {
      // std::cout << "(" << it->first << "," << it->second << "," << num_strands << ")\n";
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
  Laurent_p quantum_invariant(Plat_braid& b) {


    // std::cout << "++++ enter quantum invariant\n";

    if(b.braid().empty()) { return Laurent_p(); }

    int num_strands = b.num_strands();
    auto it = b.braid().begin();
    Morphism tau = id_R_id(it->begin(),it->end(),num_strands);
    ++it;

    // std::cout << "id R id:\n";
    // std::cout << tau << "\n\n";

    for(; it!=b.braid().end(); ++it) {
      // std::cout << "before times\n";

      // std::cout << "    layer = ";
      // for(auto pp : (*it)) {
      //   std::cout << "[" << pp.first << "," << pp.second << "]";
      // }
      // std::cout << std::endl;


      tau.ltimes_equal(id_R_id(it->begin(),it->end(),num_strands));
      // std::cout << "after times\n";

    }

    // std::cout << "before h\n";
    tau.ltimes_equal(h_tensor_);
    // std::cout << "after h\n";


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

    // std::cout << "     enter id_R_id\n";

    if(beg == end) { return id_morphism(std::pow(N,num_strands)); }

    Morphism M;//(0,0,Laurent_punction());
    auto curr_idx = beg->first;
    //initialize to \operatorname{id}_N^{\otimes (i_1-1)} \otimes R^{n_1}
    if(beg->second < 0) {//n_1 < 0 => negative twists
      M = (pow_Rinv_[std::abs(beg->second)]).ltensor(id_morphism(std::pow(N,curr_idx)));
    }
    else {//n_1 > 0 => positive twists
      M = (pow_R_[beg->second]).ltensor(id_morphism(std::pow(N,curr_idx)));
    }

    // std::cout << "     now we have: ";
    // std::cout << M << "\n";



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

      // std::cout << "     now we have: ";
      // std::cout << M << "\n";


      prev_idx = curr_idx;
      ++beg;
    }
    //add the last \otimes id_N^{num_strands-|i_k|-1}
    M.rtensor_equal(id_morphism(std::pow(N,num_strands-curr_idx-2)));

    // std::cout << "     now we have (final): ";
    // std::cout << M << "\n";


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
 * @brief      Return the quantum invariant of the unknot = x^2 - x^-2.
 *
 * @return     The quantum invariant of the unknot, as a rational function.
 */
  Laurent_p quantum_invariant_unknot() {
    return (Laurent_p(-1,1) + Laurent_p(1,1));
  }

private:
  int N;
/**
 * A certain tensor of h, fixed by the number of strands considered.
 */
  Morphism h_tensor_;

  std::vector<Morphism> pow_R_;
  std::vector<Morphism> pow_Rinv_;
};

} // namespace kumquat

#endif // END KUMQUAT_JONES_POLYNOMIAL_H_ 