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

#ifndef KUMQUAT_FREE_H_ 
#define KUMQUAT_DENSE_MATRIX_H_

#include <iomanip>
#include <vector>
#include <map>
#include <kumquat/Free_module_finite_rank.h>

namespace kumquat {

/** \brief Class implementing the standard requirement for a category of free modules (or vector spaces).
 * 
 * Objects are finite rank free modules over a ring (or finite dimensional vector spaces over a field), and matrices are morphisms are represented by matrices.
 * 
 * template CoefficientStructure is a model of Ring, and the ring on which the module is defined, and the ring of coefficients of the matrices.
 * */
template<typename CoefficientStructure>
class Category_free_modules {

  typedef CoefficientStructure Ring;
  typedef typename Ring::Element Scalar;//coefficients of matrices and vectors
  typedef Free_module_finite_rank<Ring> Object;//to satisfy requirement of Set
  typedef Object::Matrix Morphism;
  typedef Object::Vector Vector;

 /** \brief Compose a morphism to the right.
 * 
 * Set \f$\phi \rightarrow \phi \circ \psi\f$.
 * */
  void rcompose_equal(Morphism &phi, const Morphism& psi) {
    phi.rtimes_equal(psi);
  }
/** \brief Compose a morphism to the left.
 * 
 * Set \f$\phi \rightarrow \psi \circ \phi\f$.
 * */
  void lcompose_equal(Morphism &phi, const Morphism& psi) {
    phi.ltimes_equal(phi);
  }
/** \brief Compose a morphism to the right.
 * 
 * Return \f$\phi \circ \psi\f$, based on \f$\mathtt{rcompose_equal}\f$.
 * */
  Morphism compose(const Morphism &phi, const Morphism& psi) {
    return phi * psi;
  }
/** \brief Evaluate a morphism on an object.
 * 
 * Set \f$ x \leftarrow \phi(x)\f$.
 * */
  void evaluate_equal(Object &obj, const Morphism& phi);
/** \brief Evaluate a morphism on an object.
 * 
 * Return \f$ x \leftarrow \phi(x)\f$, based on \f$\mathtt{evaluate_equal}\f$.
 * */
  Object evaluate(Object &obj, const Morphism& phi);
/** \brief Test for equality of objects.*/
  bool equal(Object &o1, Object &o2) {
    return o1.rank() == o2.rank();
  }

   
};
