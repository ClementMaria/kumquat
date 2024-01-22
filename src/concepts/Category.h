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

namespace kumquat {

/** \brief Concept for a monoidal category (algebra). 
 * 
 * A MonoidalCategory defines a tensor product.
 */
struct Category {
/** \brief Type object is a model of Set.*/
  typedef unspecified Object;
/** \brief The type of arrows in the category.*/
  typedef unspecifed Morphism;
/** \brief Compose a morphism to the right.
 * 
 * Set \f$\phi \rightarrow \phi \circ \psi\f$.
 * */
  void rcompose_equal(Morphism &phi, const Morphism& psi);
/** \brief Compose a morphism to the left.
 * 
 * Set \f$\phi \rightarrow \psi \circ \phi\f$.
 * */
  void lcompose_equal(Morphism &phi, const Morphism& psi);
/** \brief Compose a morphism to the right.
 * 
 * Return \f$\phi \circ \psi\f$, based on \f$\mathtt{rcompose_equal}\f$.
 * */
  Morphism compose(const Morphism &phi, const Morphism& psi);
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
};  

} //namespace kumquat